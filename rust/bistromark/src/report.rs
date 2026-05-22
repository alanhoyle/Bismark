use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use bismark_lib::BISMARK_VERSION;

use crate::call::Outcome;

/// Alignment statistics, mirroring Bismark's %counting hash.
#[derive(Default, Debug)]
pub struct AlignStats {
    pub total_reads: u64,
    pub unique: u64,
    pub ambiguous: u64,
    pub unmapped: u64,
    pub no_genomic_seq: u64,
    pub cpg_meth: u64,
    pub cpg_unmeth: u64,
    pub chg_meth: u64,
    pub chg_unmeth: u64,
    pub chh_meth: u64,
    pub chh_unmeth: u64,
    pub unknown_meth: u64,
    pub unknown_unmeth: u64,
}

impl AlignStats {
    /// Tally XM string counts from one aligned read.
    pub fn tally_xm(&mut self, xm: &[u8]) {
        for &b in xm {
            match b {
                b'Z' => self.cpg_meth += 1,
                b'z' => self.cpg_unmeth += 1,
                b'X' => self.chg_meth += 1,
                b'x' => self.chg_unmeth += 1,
                b'H' => self.chh_meth += 1,
                b'h' => self.chh_unmeth += 1,
                b'U' => self.unknown_meth += 1,
                b'u' => self.unknown_unmeth += 1,
                _ => {}
            }
        }
    }

    pub fn tally_outcome(&mut self, outcome: Outcome) {
        self.total_reads += 1;
        match outcome {
            Outcome::Unique => self.unique += 1,
            Outcome::Ambiguous => self.ambiguous += 1,
            Outcome::Unmapped => self.unmapped += 1,
        }
    }

    fn total_c_meth(&self) -> u64 {
        self.cpg_meth + self.chg_meth + self.chh_meth
    }
    fn total_c_unmeth(&self) -> u64 {
        self.cpg_unmeth + self.chg_unmeth + self.chh_unmeth
    }

    fn pct(meth: u64, total: u64) -> String {
        if total == 0 {
            "N/A".to_owned()
        } else {
            format!("{:.1}%", meth as f64 / total as f64 * 100.0)
        }
    }
}

/// Write a SE alignment report compatible with Perl bismark's `*_SE_report.txt`.
pub fn write_se_report(
    stats: &AlignStats,
    genome_folder: &Path,
    input_path: &Path,
    report_path: &Path,
    bowtie2_cmd: &str,
) -> Result<()> {
    let mut w =
        BufWriter::new(std::fs::File::create(report_path).map_err(|e| {
            anyhow::anyhow!("cannot create report {}: {}", report_path.display(), e)
        })?);
    write_se_report_to(&mut w, stats, genome_folder, input_path, bowtie2_cmd)
}

pub fn write_se_report_to<W: Write>(
    w: &mut W,
    stats: &AlignStats,
    _genome_folder: &Path,
    input_path: &Path,
    bowtie2_cmd: &str,
) -> Result<()> {
    let total = stats.total_reads;
    let unique_pct = AlignStats::pct(stats.unique, total);

    writeln!(
        w,
        "Bismark report for: {} (version: {})",
        input_path.display(),
        BISMARK_VERSION
    )?;
    writeln!(
        w,
        "Bismark was run with Bowtie 2 and settings: {}",
        bowtie2_cmd
    )?;
    writeln!(w, "Option '--directional' specified (default mode): alignments to complementary strands will be ignored (i.e. not performed)")?;
    writeln!(w)?;
    writeln!(w, "Final Alignment report")?;
    writeln!(w, "======================")?;
    writeln!(w, "Sequences analysed in total:\t{}", total)?;
    writeln!(
        w,
        "Number of alignments with a unique best hit from the different bisulfite genomes:\t{}",
        stats.unique
    )?;
    writeln!(w, "Mapping efficiency:\t{}", unique_pct)?;
    writeln!(
        w,
        "Sequences with no alignments under any condition:\t{}",
        stats.unmapped
    )?;
    writeln!(w, "Sequences did not map uniquely:\t{}", stats.ambiguous)?;
    writeln!(
        w,
        "Sequences which were discarded because genomic sequence could not be extracted:\t{}",
        stats.no_genomic_seq
    )?;
    writeln!(w)?;
    writeln!(w, "Number of sequences with unique best (first) alignment came from the different bisulfite strands:")?;
    writeln!(w)?;
    writeln!(w, "Final Cytosine Methylation Report")?;
    writeln!(w, "=================================")?;
    let total_c = stats.total_c_meth() + stats.total_c_unmeth();
    writeln!(w, "Total number of C's analysed:\t{}", total_c)?;
    writeln!(w)?;
    writeln!(
        w,
        "Total methylated C's in CpG context:\t{}",
        stats.cpg_meth
    )?;
    writeln!(
        w,
        "Total methylated C's in CHG context:\t{}",
        stats.chg_meth
    )?;
    writeln!(
        w,
        "Total methylated C's in CHH context:\t{}",
        stats.chh_meth
    )?;
    writeln!(
        w,
        "Total methylated C's in Unknown context:\t{}",
        stats.unknown_meth
    )?;
    writeln!(w)?;
    writeln!(
        w,
        "Total unmethylated C's in CpG context:\t{}",
        stats.cpg_unmeth
    )?;
    writeln!(
        w,
        "Total unmethylated C's in CHG context:\t{}",
        stats.chg_unmeth
    )?;
    writeln!(
        w,
        "Total unmethylated C's in CHH context:\t{}",
        stats.chh_unmeth
    )?;
    writeln!(
        w,
        "Total unmethylated C's in Unknown context:\t{}",
        stats.unknown_unmeth
    )?;
    writeln!(w)?;
    writeln!(
        w,
        "C methylated in CpG context:\t{}",
        AlignStats::pct(stats.cpg_meth, stats.cpg_meth + stats.cpg_unmeth)
    )?;
    writeln!(
        w,
        "C methylated in CHG context:\t{}",
        AlignStats::pct(stats.chg_meth, stats.chg_meth + stats.chg_unmeth)
    )?;
    writeln!(
        w,
        "C methylated in CHH context:\t{}",
        AlignStats::pct(stats.chh_meth, stats.chh_meth + stats.chh_unmeth)
    )?;
    writeln!(
        w,
        "C methylated in Unknown context (CN or CHN):\t{}",
        AlignStats::pct(
            stats.unknown_meth,
            stats.unknown_meth + stats.unknown_unmeth
        )
    )?;
    Ok(())
}

/// Write a PE alignment report compatible with Perl bismark's `*_PE_report.txt`.
pub fn write_pe_report(
    stats: &AlignStats,
    genome_folder: &Path,
    r1_path: &Path,
    r2_path: &Path,
    report_path: &Path,
    bowtie2_cmd: &str,
) -> Result<()> {
    let mut w =
        BufWriter::new(std::fs::File::create(report_path).map_err(|e| {
            anyhow::anyhow!("cannot create report {}: {}", report_path.display(), e)
        })?);
    write_pe_report_to(&mut w, stats, genome_folder, r1_path, r2_path, bowtie2_cmd)
}

pub fn write_pe_report_to<W: Write>(
    w: &mut W,
    stats: &AlignStats,
    _genome_folder: &Path,
    r1_path: &Path,
    r2_path: &Path,
    bowtie2_cmd: &str,
) -> Result<()> {
    let total = stats.total_reads;
    let unique_pct = AlignStats::pct(stats.unique, total);

    writeln!(
        w,
        "Bismark report for: {} and {} (version: {})",
        r1_path.display(),
        r2_path.display(),
        BISMARK_VERSION
    )?;
    writeln!(
        w,
        "Bismark was run with Bowtie 2 and settings: {}",
        bowtie2_cmd
    )?;
    writeln!(w, "Option '--directional' specified (default mode): alignments to complementary strands will be ignored (i.e. not performed)")?;
    writeln!(w)?;
    writeln!(w, "Final Alignment report")?;
    writeln!(w, "======================")?;
    writeln!(w, "Sequence pairs analysed in total:\t{}", total)?;
    writeln!(
        w,
        "Number of paired-end alignments with a unique best hit:\t{}",
        stats.unique
    )?;
    writeln!(w, "Mapping efficiency:\t{}", unique_pct)?;
    writeln!(
        w,
        "Sequence pairs with no alignments under any condition:\t{}",
        stats.unmapped
    )?;
    writeln!(
        w,
        "Sequence pairs did not map uniquely:\t{}",
        stats.ambiguous
    )?;
    writeln!(
        w,
        "Sequence pairs which were discarded because genomic sequence could not be extracted:\t{}",
        stats.no_genomic_seq
    )?;
    writeln!(w)?;
    writeln!(w, "Number of sequence pairs with unique best (first) alignment came from the different bisulfite strands:")?;
    writeln!(w)?;
    writeln!(w, "Final Cytosine Methylation Report")?;
    writeln!(w, "=================================")?;
    let total_c = stats.total_c_meth() + stats.total_c_unmeth();
    writeln!(w, "Total number of C's analysed:\t{}", total_c)?;
    writeln!(w)?;
    writeln!(
        w,
        "Total methylated C's in CpG context:\t{}",
        stats.cpg_meth
    )?;
    writeln!(
        w,
        "Total methylated C's in CHG context:\t{}",
        stats.chg_meth
    )?;
    writeln!(
        w,
        "Total methylated C's in CHH context:\t{}",
        stats.chh_meth
    )?;
    writeln!(
        w,
        "Total methylated C's in Unknown context:\t{}",
        stats.unknown_meth
    )?;
    writeln!(w)?;
    writeln!(
        w,
        "Total unmethylated C's in CpG context:\t{}",
        stats.cpg_unmeth
    )?;
    writeln!(
        w,
        "Total unmethylated C's in CHG context:\t{}",
        stats.chg_unmeth
    )?;
    writeln!(
        w,
        "Total unmethylated C's in CHH context:\t{}",
        stats.chh_unmeth
    )?;
    writeln!(
        w,
        "Total unmethylated C's in Unknown context:\t{}",
        stats.unknown_unmeth
    )?;
    writeln!(w)?;
    writeln!(
        w,
        "C methylated in CpG context:\t{}",
        AlignStats::pct(stats.cpg_meth, stats.cpg_meth + stats.cpg_unmeth)
    )?;
    writeln!(
        w,
        "C methylated in CHG context:\t{}",
        AlignStats::pct(stats.chg_meth, stats.chg_meth + stats.chg_unmeth)
    )?;
    writeln!(
        w,
        "C methylated in CHH context:\t{}",
        AlignStats::pct(stats.chh_meth, stats.chh_meth + stats.chh_unmeth)
    )?;
    writeln!(
        w,
        "C methylated in Unknown context (CN or CHN):\t{}",
        AlignStats::pct(
            stats.unknown_meth,
            stats.unknown_meth + stats.unknown_unmeth
        )
    )?;
    Ok(())
}
