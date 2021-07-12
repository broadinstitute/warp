package cmd

import (
	"fmt"
	"os"

	"github.com/spf13/cobra"

	"github.com/spf13/viper"
)

var (
	cfgFile  string
	cacheDir string
)

const DefaultCacheDir = "~/.wreleaser/cache.json"

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "wreleaser",
	Short: "A tool to interact with WARP releases",
	Long: `wreleaser is a lightweight CLI to list, query and download various WARP releases
	
Currently available pipelines:
	- CEMBA
	- Exome Germline Sample
	- GDC Whole Genome Somatic Single Sample Array
	- Illumina Genotyping Array
	- Optimus
	- Single cell ATAC
	- Smart-seq2 Multi-Sample
	- Smart-seq2 Single Nucleus Multi-Sample
	- Smart seq2 Single Nucleus
	- Smart seq2 Single Sample
	- Whole Genome Germline Single Sample`,
}

func Execute() {
	cobra.CheckErr(rootCmd.Execute())
}

func init() {
	cobra.OnInitialize(initConfig)

	pflags := rootCmd.PersistentFlags()

	pflags.StringVar(&cfgFile, "config", "", "config file (default is $HOME/.wreleaser/config.yaml)")
	pflags.StringVar(&cacheDir, "cachedir", DefaultCacheDir, "Directory to cache release information")

	// Bind global flags
	viper.BindPFlag("cachedir", pflags.Lookup("cachedir"))
}

// initConfig reads in config file and ENV variables if set.
func initConfig() {
	if cfgFile != "" {
		// Use config file from the flag.
		viper.SetConfigFile(cfgFile)
	} else {
		// Find home directory.
		home, err := os.UserHomeDir()
		cobra.CheckErr(err)

		// Search config in home directory with name ".wreleaser" (without extension).
		viper.AddConfigPath(home)
		viper.SetConfigType("yaml")
		viper.SetConfigName(".wreleaser/config")
	}

	viper.AutomaticEnv() // read in environment variables that match

	// If a config file is found, read it in.
	if err := viper.ReadInConfig(); err == nil {
		fmt.Fprintln(os.Stderr, "Using config file:", viper.ConfigFileUsed())
	}
}
