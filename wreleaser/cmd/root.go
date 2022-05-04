package cmd

import (
	"os"
	"path/filepath"

	"github.com/spf13/cobra"

	"github.com/spf13/viper"
)

var (
	cacheDir   string
	outputFile string
)

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "wreleaser",
	Short: "A tool to interact with WARP releases",
	Long: `wreleaser is a lightweight CLI to list, query and download various WARP releases

Currently available pipelines:
	- AnnotationFiltration
	- Arrays
	- CEMBA
	- CramToUnmappedBams
	- ExomeGermlineSingleSample
	- ExomeReprocessing
	- ExternalExomeReprocessing
	- ExternalWholeGenomeReprocessing
	- GDCWholeGenomeSomaticSingleSample
	- IlluminaGenotypingArray
        - Imputation
	- JointGenotyping
        - JointGenotypingByChromosomePartOne
        - JointGenotypingByChromosomePartTwo
        - MultiSampleArrays
        - MultiSampleSmartSeq2
        - MultiSampleSmartSeq2SingleNucleus
        - Optimus
        - ReblockGVCF
        - RNAWithUMIsPipeline
        - SmartSeq2SingleNucleus
        - SmartSeq2SingleSample
        - ValidateChip
        - VariantCalling
        - WholeGenomeGermlineSingleSample
        - WholeGenomeReprocessing`,
}

// Execute is the main function
func Execute() {
	cobra.CheckErr(rootCmd.Execute())
}

func init() {
	pflags := rootCmd.PersistentFlags()

	pflags.StringVar(&cacheDir, "cachedir", defaultCacheDir(), "Directory to cache release information")
	pflags.StringVar(&outputFile, "output", "", "Pipe command output to file")

	// Bind global flags
	viper.BindPFlag("cachedir", pflags.Lookup("cachedir"))
	viper.BindPFlag("output", pflags.Lookup("output"))
}

func defaultCacheDir() string {
	home, err := os.UserHomeDir()
	if err != nil {
		panic("Error locating user's $HOME")
	}
	return filepath.Join(home, "/.wreleaser/cache.json")
}
