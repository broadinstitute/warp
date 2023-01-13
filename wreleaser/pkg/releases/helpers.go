package releases

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/spf13/viper"
)

var validReleases = []string{
	"AnnotationFiltration",
	"Arrays",
	"BroadInternalRNAWithUMIs",
	"BroadInternalArrays",
	"BroadInternalImputation",
	"BroadInternalUltimaGenomics",
	"CEMBA",
	"CheckFingerprint",
	"CramToUnmappedBams",
	"ExomeGermlineSingleSample",
	"ExomeReprocessing",
	"ExternalExomeReprocessing",
	"ExternalWholeGenomeReprocessing",
	"GDCWholeGenomeSomaticSingleSample",
	"IlluminaGenotypingArray",
	"Imputation",
	"JointGenotyping",
	"JointGenotypingByChromosomePartOne",
	"JointGenotypingByChromosomePartTwo",
	"MultiSampleArrays",
	"MultiSampleSmartSeq2",
	"MultiSampleSmartSeq2SingleNucleus",
	"Optimus",
	"ReblockGVCF",
	"RNAWithUMIsPipeline",
	"SmartSeq2SingleNucleus",
	"SmartSeq2SingleSample",
	"UltimaGenomicsWholeGenomeGermline",
	"UltimaGenomicsJointGenotyping",
	"ValidateChip",
	"VariantCalling",
	"WholeGenomeGermlineSingleSample",
	"WholeGenomeReprocessing",
}

var isValidPipeline = false

func cacheExists() bool {
	cacheFile := getCacheDir()

	info, err := os.Stat(cacheFile)
	if os.IsNotExist(err) {
		return false
	}
	return !info.IsDir()
}

func makeCache() (*os.File, error) {
	cacheFile := getCacheDir()

	if err := os.MkdirAll(filepath.Dir(cacheFile), 0770); err != nil {
		return nil, err
	}
	return os.Create(cacheFile)
}

func getCacheDir() string {
	return viper.GetString("cachedir")
}

// hasValue checks to see if the current pipeline exists in requested arguments
func hasValue(requested []string, pipeline string) bool {
	for _, element := range requested {
		if element == pipeline {
			return true
		}
	}
	return false
}

// ValidateArgs verifys that the requested pipelines have been released in WARP
func ValidateArgs(requested []string) {
	for _, r := range requested {
		for _, v := range validReleases {
			if r == v {
				isValidPipeline = true
			}
		}
		if !isValidPipeline {
			fmt.Fprintf(os.Stderr, "Requested pipeline %s does not exist in WARP releases \n", r)
			fmt.Fprintf(os.Stderr, "Run 'wreleaser --help' for a list of valid pipelines \n")
			os.Exit(1)
		}
		isValidPipeline = false
	}
	return
}
