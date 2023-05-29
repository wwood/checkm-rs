extern crate csv;
extern crate env_logger;
extern crate log;

use log::*;

use std::collections::BTreeMap;

pub struct CheckM1TabTable {}

impl CheckM1TabTable {
    pub fn good_quality_genome_names(
        file_path: &str,
        min_completeness: f32,
        max_contamination: f32,
    ) -> Result<Vec<String>, CheckMReadError> {
        let mut passes = vec![];
        let qualities = CheckM1TabTable::read_file_path(file_path)?;
        for (genome, quality) in qualities.genome_to_quality.iter() {
            trace!("Genome: {}, Quality: {:?}", genome, quality);
            if quality.completeness >= min_completeness
                && quality.contamination <= max_contamination
            {
                passes.push(genome.clone())
            }
        }
        debug!(
            "Read in {} genomes from {}, {} passed the quality thresholds",
            qualities.genome_to_quality.len(),
            file_path,
            passes.len()
        );
        Ok(passes)
    }

    pub fn read_file_path(
        file_path: &str,
    ) -> Result<CheckMResult<CheckM1GenomeQuality>, CheckMReadError> {
        let mut qualities = BTreeMap::new();
        let rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(std::path::Path::new(file_path));
        let mut total_seen = 0usize;

        if rdr.is_err() {
            return Err(CheckMReadError {
                msg: format!("Failed to parse CheckM v1 genome quality {}", file_path),
            });
        }

        for result in rdr.unwrap().records() {
            let res = result.expect("Parsing error in CheckM tab table file");
            if res.len() != 14 {
                return Err(CheckMReadError {
                    msg: format!(
                        "Parsing error in CheckM tab table file - didn't find 14 columns in line {:?}",
                        res
                    ),
                });
            }
            let completeness: f32 = res[11]
                .parse::<f32>()
                .expect("Error parsing completeness in checkm tab table");
            let contamination: f32 = res[12]
                .parse::<f32>()
                .expect("Error parsing contamination in checkm tab table");
            let strain_heterogeneity: f32 = res[13]
                .parse::<f32>()
                .expect("Error parsing contamination in checkm tab table");
            trace!(
                "For {}, found completeness {} and contamination {}",
                &res[0],
                completeness,
                contamination
            );
            match qualities.insert(
                res[0].to_string(),
                CheckM1GenomeQuality {
                    completeness: completeness / 100.,
                    contamination: contamination / 100.,
                    strain_heterogeneity: strain_heterogeneity / 100.,
                },
            ) {
                None => {}
                Some(_) => {
                    return Err(CheckMReadError {
                        msg: format!(
                            "The genome {} was found multiple times in the checkm file {}",
                            &res[0], file_path
                        ),
                    });
                }
            };
            total_seen += 1;
        }
        debug!("Read in {} genomes from {}", total_seen, file_path);
        Ok(CheckMResult {
            genome_to_quality: qualities,
        })
    }
}
pub struct CheckM2QualityReport {}
impl CheckM2QualityReport {
    pub fn good_quality_genome_names(
        file_path: &str,
        min_completeness: f32,
        max_contamination: f32,
    ) -> Result<Vec<String>, CheckMReadError> {
        let mut passes = vec![];
        let qualities = CheckM2QualityReport::read_file_path(file_path)?;
        for (genome, quality) in qualities.genome_to_quality.iter() {
            trace!("Genome: {}, Quality: {:?}", genome, quality);
            if quality.completeness >= min_completeness
                && quality.contamination <= max_contamination
            {
                passes.push(genome.clone())
            }
        }
        debug!(
            "Read in {} genomes from {}, {} passed the quality thresholds",
            qualities.genome_to_quality.len(),
            file_path,
            passes.len()
        );
        Ok(passes)
    }

    pub fn read_file_path(
        file_path: &str,
    ) -> Result<CheckMResult<CheckM2GenomeQuality>, CheckMReadError> {
        let mut qualities = BTreeMap::new();
        let rdr_res = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(std::path::Path::new(file_path));
        let mut total_seen = 0usize;

        if rdr_res.is_err() {
            return Err(CheckMReadError {
                msg: format!("Failed to parse CheckM v2 genome quality {}", file_path),
            });
        }
        let mut rdr = rdr_res.unwrap();

        // The number of columns can change, so just check the first 3 columns have the expected names
        let headers = rdr.headers().unwrap().clone();
        if &headers[0] != "Name" || &headers[1] != "Completeness" || &headers[2] != "Contamination"
        {
            return Err(CheckMReadError {
                msg: format!(
                    "Parsing error in CheckM2 qualities file - didn't find expected headers in line {:?}",
                    headers
                ),
            });
        }

        for result in rdr.records() {
            let res = result.expect("Parsing error in CheckM qualities file");
            let completeness: f32 = res[1]
                .parse::<f32>()
                .expect("Error parsing completeness in checkm tab table");
            let contamination: f32 = res[2]
                .parse::<f32>()
                .expect("Error parsing contamination in checkm tab table");
            trace!(
                "For {}, found completeness {} and contamination {}",
                &res[0],
                completeness,
                contamination
            );
            match qualities.insert(
                res[0].to_string(),
                CheckM2GenomeQuality {
                    completeness: completeness / 100.,
                    contamination: contamination / 100.,
                },
            ) {
                None => {}
                Some(_) => {
                    return Err(CheckMReadError {
                        msg: format!(
                            "The genome {} was found multiple times in the checkm file {}",
                            &res[0], file_path
                        ),
                    });
                }
            };
            total_seen += 1;
        }
        debug!("Read in {} genomes from {}", total_seen, file_path);
        Ok(CheckMResult {
            genome_to_quality: qualities,
        })
    }
}

pub trait GenomeQuality {
    fn completeness(&self) -> f32;
    fn contamination(&self) -> f32;
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct CheckM1GenomeQuality {
    pub completeness: f32,
    pub contamination: f32,
    pub strain_heterogeneity: f32,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct CheckM2GenomeQuality {
    pub completeness: f32,
    pub contamination: f32,
}

pub struct CheckMResult<T: GenomeQuality> {
    pub genome_to_quality: BTreeMap<String, T>,
}

impl GenomeQuality for CheckM1GenomeQuality {
    fn completeness(&self) -> f32 {
        self.completeness
    }
    fn contamination(&self) -> f32 {
        self.contamination
    }
}

impl GenomeQuality for CheckM2GenomeQuality {
    fn completeness(&self) -> f32 {
        self.completeness
    }
    fn contamination(&self) -> f32 {
        self.contamination
    }
}

impl<T: GenomeQuality + Copy + std::fmt::Debug> CheckMResult<T> {
    pub fn order_genomes_by_completeness_minus_4contamination(&self) -> Vec<&String> {
        let mut genomes_and_qualities: Vec<_> = self
            .genome_to_quality
            .iter()
            .map(|(genome, quality)| (genome, quality))
            .collect();

        // sort in reverse order
        genomes_and_qualities.sort_unstable_by(|(_, q1), (_, q2)| {
            (q2.completeness() - 4. * q2.contamination())
                .partial_cmp(&(q1.completeness() - 4. * q1.contamination()))
                .unwrap()
        });
        genomes_and_qualities
            .into_iter()
            .map(|(genome, _)| genome)
            .collect()
    }

    pub fn order_genomes_by_completeness_minus_5contamination(&self) -> Vec<&String> {
        let mut genomes_and_qualities: Vec<_> = self
            .genome_to_quality
            .iter()
            .map(|(genome, quality)| (genome, quality))
            .collect();

        // sort in reverse order
        genomes_and_qualities.sort_unstable_by(|(_, q1), (_, q2)| {
            (q2.completeness() - 5. * q2.contamination())
                .partial_cmp(&(q1.completeness() - 5. * q1.contamination()))
                .unwrap()
        });
        genomes_and_qualities
            .into_iter()
            .map(|(genome, _)| genome)
            .collect()
    }

    /// Map paths to FASTA paths to CheckM qualities, and return paths ordered
    /// by their quality, where quality is completeness - 4*contamination. If
    /// not None, min_completeness and max_contamination specify thresholds as
    /// fractions e.g. 0.8 not 80.
    pub fn order_fasta_paths_by_completeness_minus_4contamination<'a>(
        &self,
        genome_fasta_files: &[&'a str],
        min_completeness: Option<f32>,
        max_contamination: Option<f32>,
    ) -> Result<Vec<&'a str>, String> {
        let mut key_to_order = BTreeMap::new();
        for (i, key) in self
            .order_genomes_by_completeness_minus_4contamination()
            .into_iter()
            .enumerate()
        {
            key_to_order.insert(key.as_str(), i);
        }
        self.order_fasta_list(
            &key_to_order,
            genome_fasta_files,
            min_completeness,
            max_contamination,
        )
    }

    /// Map paths to FASTA paths to CheckM qualities, and return paths ordered
    /// by their quality, where quality is completeness - 5*contamination. If
    /// not None, min_completeness and max_contamination specify thresholds as
    /// fractions e.g. 0.8 not 80.
    pub fn order_fasta_paths_by_completeness_minus_5contamination<'a>(
        &self,
        genome_fasta_files: &[&'a str],
        min_completeness: Option<f32>,
        max_contamination: Option<f32>,
    ) -> Result<Vec<&'a str>, String> {
        let mut key_to_order = BTreeMap::new();
        for (i, key) in self
            .order_genomes_by_completeness_minus_5contamination()
            .into_iter()
            .enumerate()
        {
            key_to_order.insert(key.as_str(), i);
        }
        self.order_fasta_list(
            &key_to_order,
            genome_fasta_files,
            min_completeness,
            max_contamination,
        )
    }

    fn order_fasta_list<'a>(
        &self,
        key_to_order: &BTreeMap<&str, usize>,
        genome_fasta_files: &[&'a str],
        min_completeness: Option<f32>,
        max_contamination: Option<f32>,
    ) -> Result<Vec<&'a str>, String> {
        let mut fasta_and_order: Vec<(&str, usize)> = vec![];
        for fasta_path in genome_fasta_files.iter() {
            let checkm_name = std::path::Path::new(fasta_path)
                .file_stem()
                .unwrap()
                .to_str()
                .unwrap();
            match key_to_order.get(checkm_name) {
                Some(rank) => {
                    if (min_completeness.is_none()
                        || self.genome_to_quality[checkm_name].completeness()
                            >= min_completeness.unwrap())
                        && (max_contamination.is_none()
                            || self.genome_to_quality[checkm_name].contamination()
                                <= max_contamination.unwrap())
                    {
                        fasta_and_order.push((*fasta_path, *rank))
                    }
                }
                None => {
                    return Err(format!(
                        "Unable to find quality for genome fasta file {}",
                        fasta_path
                    ));
                }
            }
        }
        fasta_and_order.sort_unstable_by(|a, b| a.1.cmp(&b.1));
        Ok(fasta_and_order.into_iter().map(|(g, _)| g).collect())
    }

    pub fn filter(&self, min_completeness: f32, max_contamination: f32) -> CheckMResult<T> {
        let mut new = BTreeMap::new();

        for (g, q) in self.genome_to_quality.iter() {
            if q.completeness() >= min_completeness && q.contamination() <= max_contamination {
                new.insert(g.clone().to_string(), *q);
            }
        }
        CheckMResult {
            genome_to_quality: new,
        }
    }

    pub fn retrieve_via_fasta_path(&self, fasta_path: &str) -> Result<T, CheckMRetrievalError> {
        let checkm_name_stem_res = std::path::Path::new(fasta_path).file_stem();
        if checkm_name_stem_res.is_none() {
            return Err(CheckMRetrievalError {
                msg: format!("Unable to find file_stem for {}", fasta_path),
            });
        }
        let checkm_name_stem = checkm_name_stem_res.unwrap();
        let checkm_name_res = checkm_name_stem.to_str();
        if checkm_name_res.is_none() {
            return Err(CheckMRetrievalError {
                msg: format!(
                    "Failed to convert fasta file name to string: {}",
                    fasta_path
                ),
            });
        }
        let checkm_name = checkm_name_res.unwrap();
        debug!(
            "Retrieving checkm name {}, derived from {}",
            checkm_name, fasta_path
        );
        trace!("{:?}", self.genome_to_quality);
        match self.genome_to_quality.get(checkm_name) {
            Some(q) => Ok(*q),
            None => {
                // Possibly the checkm file was created with absolute paths. Try
                // that.
                let checkm_parent_res = std::path::Path::new(fasta_path).parent();
                if checkm_parent_res.is_none() {
                    return Err(CheckMRetrievalError {
                        msg: format!("Unable to find parent of fasta path {}", fasta_path),
                    });
                };
                let joined = checkm_parent_res.unwrap().join(checkm_name_stem);
                let checkm_name2 = joined.to_str();
                if checkm_name2.is_none() {
                    return Err(CheckMRetrievalError {
                        msg: format!("Unable to convert name {} to str", checkm_name),
                    });
                }
                debug!(
                    "Retrieving absolute path checkm name {}, derived from {}",
                    checkm_name2.unwrap(),
                    fasta_path
                );
                match self.genome_to_quality.get(checkm_name2.unwrap()) {
                    Some(q) => Ok(*q),
                    None => Err(CheckMRetrievalError {
                        msg: format!(
                            "Unable to find checkm name {} in checkm results",
                            checkm_name
                        ),
                    }),
                }
            }
        }
    }
}
#[derive(Debug, PartialEq)]
pub struct CheckMRetrievalError {
    msg: String,
}

impl std::fmt::Display for CheckMRetrievalError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.msg)
    }
}

impl std::error::Error for CheckMRetrievalError {}

#[derive(Debug, PartialEq)]
pub struct CheckMReadError {
    msg: String,
}

impl std::fmt::Display for CheckMReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.msg)
    }
}

impl std::error::Error for CheckMReadError {}

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
            vec![
                "GUT_GENOME006390.gff",
                "GUT_GENOME011264.gff",
                "GUT_GENOME011296.gff",
                "GUT_GENOME011536.gff"
            ],
            CheckM1TabTable::good_quality_genome_names(&"tests/data/checkm.tsv", 0.56, 0.1)
                .unwrap()
        )
    }

    #[test]
    fn test_checkm2_good_quality_genome_names() {
        init();
        assert_eq!(
            Ok(vec!["UTC-1_IN_HT_bin.22".to_string(),]),
            CheckM2QualityReport::good_quality_genome_names(
                &"tests/data/checkm2/quality_report.tsv",
                0.7,
                0.1
            )
        );
        let empty: Vec<String> = vec![];
        assert_eq!(
            Ok(empty),
            CheckM2QualityReport::good_quality_genome_names(
                &"tests/data/checkm2/quality_report.tsv",
                0.75,
                0.1
            )
        )
    }

    #[test]
    fn test_checkm2_14column_format() {
        init();
        assert_eq!(
            Ok(vec!["AAMD-1".to_string(),]),
            CheckM2QualityReport::good_quality_genome_names(
                &"tests/data/checkm2/quality_report2.tsv",
                0.7,
                0.1
            )
        );
    }

    #[test]
    fn test_retrieve() {
        init();
        let checkm = CheckM1TabTable::read_file_path(&"tests/data/checkm.tsv").unwrap();
        assert_eq!(
            Ok(CheckM1GenomeQuality {
                completeness: 83.38 / 100.,
                contamination: 0.,
                strain_heterogeneity: 0.
            }),
            checkm.retrieve_via_fasta_path(&"/some/path/GUT_GENOME011264.gff.fna")
        );
        assert!(checkm
            .retrieve_via_fasta_path(&"/some/path/GUT_not_a_genome_GENOME011264.gff.fna")
            .is_err());
    }

    #[test]
    fn test_ordering_4times() {
        init();
        let checkm = CheckM1TabTable::read_file_path(&"tests/data/checkm.tsv").unwrap();
        assert_eq!(
            vec![
                "GUT_GENOME006390.gff",
                "GUT_GENOME011264.gff",
                "GUT_GENOME011296.gff",
                "GUT_GENOME011367.gff",
                "GUT_GENOME011536.gff"
            ],
            checkm.order_genomes_by_completeness_minus_4contamination()
        );
        let checkm = CheckM1TabTable::read_file_path(&"tests/data/checkm2.tsv").unwrap();
        assert_eq!(
            vec![
                "GUT_GENOME006390.gff",
                "GUT_GENOME011264.gff",
                "GUT_GENOME011296.gff",
                "GUT_GENOME011536.gff",
                "GUT_GENOME011367.gff"
            ],
            checkm.order_genomes_by_completeness_minus_4contamination()
        );
    }

    #[test]
    fn test_ordering_5times() {
        init();
        // Cheating here a bit - same result as the minus_4times
        let checkm = CheckM1TabTable::read_file_path(&"tests/data/checkm.tsv").unwrap();
        assert_eq!(
            vec![
                "GUT_GENOME006390.gff",
                "GUT_GENOME011264.gff",
                "GUT_GENOME011296.gff",
                "GUT_GENOME011367.gff",
                "GUT_GENOME011536.gff"
            ],
            checkm.order_genomes_by_completeness_minus_5contamination()
        );
        let checkm = CheckM1TabTable::read_file_path(&"tests/data/checkm2.tsv").unwrap();
        assert_eq!(
            vec![
                "GUT_GENOME006390.gff",
                "GUT_GENOME011264.gff",
                "GUT_GENOME011296.gff",
                "GUT_GENOME011536.gff",
                "GUT_GENOME011367.gff"
            ],
            checkm.order_genomes_by_completeness_minus_5contamination()
        );
    }

    #[test]
    fn test_fasta_ordering_4times() {
        init();
        let checkm = CheckM1TabTable::read_file_path(&"tests/data/checkm.tsv").unwrap();
        assert_eq!(
            vec![
                "/tmp/GUT_GENOME006390.gff.fna",
                "GUT_GENOME011264.gff.fna",
                "GUT_GENOME011296.gff.fna",
                "GUT_GENOME011367.gff.fna",
                "GUT_GENOME011536.gff.fna"
            ],
            checkm
                .order_fasta_paths_by_completeness_minus_4contamination(
                    &vec![
                        "GUT_GENOME011264.gff.fna",
                        "/tmp/GUT_GENOME006390.gff.fna",
                        "GUT_GENOME011296.gff.fna",
                        "GUT_GENOME011367.gff.fna",
                        "GUT_GENOME011536.gff.fna"
                    ],
                    None,
                    None,
                )
                .unwrap()
        );
        let checkm = CheckM1TabTable::read_file_path(&"tests/data/checkm2.tsv").unwrap();
        assert_eq!(
            vec![
                "/tmp/GUT_GENOME006390.gff.fna",
                "GUT_GENOME011264.gff.fna",
                "GUT_GENOME011296.gff.fna",
                "GUT_GENOME011536.gff.fna",
                "GUT_GENOME011367.gff.fna"
            ],
            checkm
                .order_fasta_paths_by_completeness_minus_4contamination(
                    &vec![
                        "GUT_GENOME011264.gff.fna",
                        "/tmp/GUT_GENOME006390.gff.fna",
                        "GUT_GENOME011296.gff.fna",
                        "GUT_GENOME011367.gff.fna",
                        "GUT_GENOME011536.gff.fna"
                    ],
                    None,
                    None,
                )
                .unwrap()
        );
    }

    #[test]
    fn test_fasta_ordering_4times_min_completeness() {
        init();
        let checkm = CheckM1TabTable::read_file_path(&"tests/data/checkm.tsv").unwrap();
        assert_eq!(
            vec![
                "/tmp/GUT_GENOME006390.gff.fna",
                "GUT_GENOME011264.gff.fna",
                "GUT_GENOME011296.gff.fna"
            ],
            checkm
                .order_fasta_paths_by_completeness_minus_4contamination(
                    &vec![
                        "GUT_GENOME011264.gff.fna",
                        "/tmp/GUT_GENOME006390.gff.fna",
                        "GUT_GENOME011296.gff.fna",
                        "GUT_GENOME011367.gff.fna",
                        "GUT_GENOME011536.gff.fna"
                    ],
                    Some(0.7),
                    None,
                )
                .unwrap()
        );
    }

    #[test]
    fn test_absolute_path_retrieval() {
        init();
        let checkm = CheckM1TabTable::read_file_path(&"tests/data/checkm3.tsv").unwrap();
        assert_eq!(
            Ok(CheckM1GenomeQuality {
                completeness: 0.9361,
                contamination: 0.0037,
                strain_heterogeneity: 1.0,
            }),
            checkm.retrieve_via_fasta_path("/tmp/GUT_GENOME006390.gff.fna")
        );
    }
}
