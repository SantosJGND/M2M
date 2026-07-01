use clap::{ArgAction, Parser};

#[derive(Parser, Debug)]
#[command(
    version,
    about,
    long_about = "a tool to compare mappings from multiple files"
)]
pub struct Cli {
    /// files to parse, comma separated
    #[clap(long, default_value = "")]
    pub files: String,
    pub _files_list: Vec<String>,

    /// input directory
    #[clap(short, long, default_value = "")]
    pub input_directory: String,

    /// maximum reads to keep
    #[clap(short, long, default_value = "500000")]
    pub max_reads: i32,

    // keep samples with no reads mapped in distance matrix (will be filtered out in tree building).
    // Can give counter-intuitive results in symmetric cases.
    #[clap(long, action = ArgAction::SetTrue, default_value_t = false)]
    pub _keep_empty_samples: bool,

    /// use adjusted distances in neighbour joining algorithm.
    /// Can give counter-intuitive results in symmetric cases.
    #[clap(short, long, default_value = "false")]
    pub adjusted_distances: bool,

    /// private reads proportion threshold.
    /// Nodes with private proportion > threshold will be reported.
    /// Default: 0.6
    #[clap(short, long, default_value = "0.6")]
    pub threshold: f64,

    /// minimum frequency threshold - calculated across files.
    #[clap(short, long, default_value = "0.1")]
    pub frequency_threshold: f32,

    /// output directory    
    /// default: output
    #[clap(short, long, default_value = "output")]
    pub output_directory: String,

    /// enable cluster analysis (tree building and clade identification).
    /// When disabled, only distance matrices are generated.
    /// Use --no-cluster-analysis to disable.
    #[clap(long, action = ArgAction::SetTrue, default_value_t = true)]
    #[clap(long = "no-cluster-analysis", action = ArgAction::SetFalse)]
    pub cluster_analysis: bool,
}

impl Cli {
    pub fn new() -> Self {
        let args = Cli::parse();
        Cli::parse_args(args)
    }

    fn parse_args(args: Cli) -> Cli {
        let files = args.files;

        let _files_list = files
            .split(",")
            .map(|s| s.to_string())
            .filter(|s| !s.is_empty())
            .collect();

        Cli {
            files,
            _files_list,
            input_directory: args.input_directory,
            max_reads: args.max_reads,
            _keep_empty_samples: args._keep_empty_samples,
            adjusted_distances: args.adjusted_distances,
            threshold: args.threshold,
            frequency_threshold: args.frequency_threshold,
            output_directory: args.output_directory,
            cluster_analysis: args.cluster_analysis,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cli() {
        let cli = Cli::new();
        assert_eq!(cli.max_reads, 500000);
        assert_eq!(cli.adjusted_distances, false);
        assert_eq!(cli.frequency_threshold, 0.1);
        assert_eq!(cli.threshold, 0.6);
        assert_eq!(cli.output_directory, "output");
        assert_eq!(cli.cluster_analysis, true);
    }
}
