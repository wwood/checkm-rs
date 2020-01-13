extern crate log;
extern crate env_logger;
extern crate csv;

use std;
use std::collections::BTreeMap;
use log::*;

pub struct CheckMTabTable {
}

impl CheckMTabTable {
    pub fn good_quality_genome_names(file_path: &str, min_completeness: f32, max_contamination: f32) -> Vec<String> {
        let mut passes = vec![];
        let qualities = CheckMTabTable::read_file_path(file_path);
        for (genome, quality) in qualities.genome_to_quality.iter() {
            trace!("Genome: {}, Quality: {:?}", genome, quality);
            if quality.completeness >= min_completeness && quality.contamination <= max_contamination {
                passes.push(genome.clone())
            }
        }
        debug!("Read in {} genomes from {}, {} passed the quality thresholds", 
            qualities.genome_to_quality.len(), file_path, passes.len());
        return passes;
    }

    pub fn read_file_path(file_path: &str) -> CheckMResult {
        let mut qualities = BTreeMap::new();
        let rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(std::path::Path::new(file_path));
        let mut total_seen = 0usize;

        for result in rdr
            .expect(&format!("Failed to parse CheckM tab table {}", file_path))
            .records() {

            let res = result.expect("Parsing error in CheckM tab table file");
            if res.len() != 14 {
                error!("Parsing error in CheckM tab table file - didn't find 13 columns in line {:?}", res);
                std::process::exit(1);
            }
            let completeness: f32 = res[11].parse::<f32>().expect("Error parsing completeness in checkm tab table");
            let contamination: f32 = res[12].parse::<f32>().expect("Error parsing contamination in checkm tab table");
            trace!("For {}, found completeness {} and contamination {}", &res[0], completeness, contamination);
            match qualities.insert(
                res[0].to_string(),
                GenomeQuality { 
                    completeness: completeness/100.,
                    contamination: contamination/100.,
                }) {
                None => {},
                Some(_) => {
                    error!("The genome {} was found multiple times in the checkm file {}", res[0].to_string(), file_path);
                    std::process::exit(1);
                }
            };
            total_seen += 1;
        }
        debug!("Read in {} genomes from {}", total_seen, file_path);
        return CheckMResult {
            genome_to_quality: qualities
        };
    }
}

pub struct CheckMResult {
    pub genome_to_quality: BTreeMap<String,GenomeQuality>,
}

#[derive(Debug,PartialEq)]
pub struct GenomeQuality {
    pub completeness: f32,
    pub contamination: f32,
}

impl CheckMResult {
    pub fn retrieve_via_fasta_path(&self, fasta_path: &str) -> Result<&GenomeQuality,()> {
        let checkm_name = std::path::Path::new(fasta_path).file_stem().unwrap().to_str().unwrap();
        match self.genome_to_quality.get(checkm_name) {
            Some(q) => Ok(q),
            None => Err(())
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }
    
    #[test]
    fn test_good_quality_genome_names() {
        init();
        assert_eq!(
            vec!["GUT_GENOME006390.gff","GUT_GENOME011264.gff","GUT_GENOME011296.gff","GUT_GENOME011536.gff"], 
            CheckMTabTable::good_quality_genome_names(&"tests/data/checkm.tsv",0.56,0.1))
    }

    #[test]
    fn test_retrieve() {
        init();
        let checkm = CheckMTabTable::read_file_path(&"tests/data/checkm.tsv");
        assert_eq!(
            Ok(&GenomeQuality {
                completeness: 83.38/100.,
                contamination: 0.,
            }),
            checkm.retrieve_via_fasta_path(&"/some/path/GUT_GENOME011264.gff.fna"));
        assert_eq!(
            Err(()),
            checkm.retrieve_via_fasta_path(&"/some/path/GUT_not_a_genome_GENOME011264.gff.fna"));
    }
}