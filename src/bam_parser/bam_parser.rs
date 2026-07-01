use bam::header::Header;
use bam::record::Record;
use bio::io::fastq;
use std::collections::HashSet;
use std::io;

pub struct BamParser {
    pub filepath: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct HashedFastqRecord {
    pub id: String,
    pub seq: String,
    pub qual: String,
}

impl HashedFastqRecord {
    pub fn from_record(record: &fastq::Record) -> Self {
        HashedFastqRecord {
            id: record.id().to_string(),
            seq: String::from_utf8_lossy(record.seq()).to_string(),
            qual: String::from_utf8_lossy(record.qual()).to_string(),
        }
    }

    pub fn new(id: String, seq: String, qual: String) -> Self {
        HashedFastqRecord { id, seq, qual }
    }
}

impl BamParser {
    pub fn new(filepath: String) -> Self {
        BamParser { filepath }
    }

    pub fn parse_header(&self) -> io::Result<Header> {
        let reader = bam::BamReader::from_path(&self.filepath, 4)?;
        let header = reader.header().clone();
        Ok(header)
    }

    pub fn parse_records(&self) -> io::Result<Vec<Record>> {
        let mut reader = bam::BamReader::from_path(&self.filepath, 4)?;
        let mut records = Vec::new();
        while let Some(result) = reader.next() {
            let record = result?;
            records.push(record);
        }
        Ok(records)
    }

    pub fn collect_fastq_records(&self) -> io::Result<Vec<fastq::Record>> {
        let mut fastq_hash_records = HashSet::new();
        let mut reader = bam::BamReader::from_path(&self.filepath, 4)?;
        while let Some(result) = reader.next() {
            let bam_record = result?;

            // ensure that the sequence is in the correct format
            let fastq_record = fastq::Record::with_attrs(
                String::from_utf8_lossy(bam_record.name())
                    .to_string()
                    .as_str(),
                Some(""),
                &bam_record.sequence().to_vec(),
                &bam_record.qualities().to_readable(),
            );

            let hashed_record = HashedFastqRecord::from_record(&fastq_record);

            if !fastq_hash_records.contains(&hashed_record) {
                fastq_hash_records.insert(hashed_record);
            }
        }
        let mut fastq_records = Vec::new();
        for record in fastq_hash_records {
            let fastq_record = fastq::Record::with_attrs(
                record.id.as_str(),
                Some("".to_string().as_str()),
                &record.seq.into_bytes(),
                &record.qual.into_bytes(),
            );
            fastq_records.push(fastq_record);
        }

        Ok(fastq_records)
    }

    pub fn write_fastq_records(
        &self,
        fastq_records: Vec<fastq::Record>,
        output_path: &str,
    ) -> io::Result<()> {
        let mut writer = fastq::Writer::to_file(output_path)?;
        for record in fastq_records {
            writer.write_record(&record)?;
        }
        writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::collections::HashSet;

    #[test]
    fn test_bam_parser() {
        let bam_file = "/home/bioinf/Desktop/CODE/PERSONAL/TOOlS/bam_parser/data/vg230245_rpip_run5.JX627336_1.sorted.bam"; // Replace with a valid BAM file path
        let parser = BamParser::new(bam_file.to_string());
        let records = parser.parse_records().expect("Failed to parse records");
        assert!(!records.is_empty());
    }

    #[test]
    fn test_collect_fastq_records() {
        let bam_file = "/home/bioinf/Desktop/CODE/PERSONAL/TOOlS/bam_parser/data/vg230245_rpip_run5.JX627336_1.sorted.bam"; // Replace with a valid BAM file path
        let parser = BamParser::new(bam_file.to_string());
        let fastq_records = parser
            .collect_fastq_records()
            .expect("Failed to collect FASTQ records");
        assert!(!fastq_records.is_empty());
        let unique_records: HashSet<_> = fastq_records.iter().collect();
        assert_eq!(
            fastq_records.len(),
            unique_records.len(),
            "There are duplicate FASTQ records"
        );
    }

    #[test]
    fn test_write_fastq_records() {
        let bam_file = "/home/bioinf/Desktop/CODE/PERSONAL/TOOlS/bam_parser/data/vg230245_rpip_run5.JX627336_1.sorted.bam"; // Replace with a valid BAM file path
        let parser = BamParser::new(bam_file.to_string());
        let fastq_records = parser
            .collect_fastq_records()
            .expect("Failed to collect FASTQ records");
        let output_path: &'static str = "output.fastq";
        parser
            .write_fastq_records(fastq_records, output_path)
            .expect("Failed to write FASTQ records");

        // Verify that the output file was created
        let output_file = std::fs::File::open(output_path);
        assert!(output_file.is_ok(), "Output file was not created");
        // Clean up the output file after the test
        std::fs::remove_file(output_path).expect("Failed to remove output file");
    }
}
