package cmd

import (
	"github.com/spf13/cobra"

	"github.com/spf13/viper"

	"github.com/broadinstitute/warp/wreleaser/pkg/releases"
)

var verbosity bool

// infoCmd represents the info command
var infoCmd = &cobra.Command{
	Use:   "info",
	Short: "Display information for specified pipeline(s)",
	Long: `This command will display information for specified pipeline(s)

Default information (non-verbose) provided for the pipeline:
	- Release version
	- Release notes
	- Release URL

Usage examples:

'wreleaser info' (display default information for all available pipelines)

'wreleaser info -v' (display versbose information for all available pipelines)

'wreleaser info CEMBA ExomeGermlineSample Optimus' (display info for the CEMBA, ExomeGermline and Optimus pipelines)

`,
	Run: func(cmd *cobra.Command, args []string) {
		releases.NewReleaseList()
	},
}

func init() {
	rootCmd.AddCommand(infoCmd)

	pflags := infoCmd.PersistentFlags()

	pflags.BoolVarP(&verbosity, "verbosity", "v", false, "Logging verbosity (default false)")

	viper.BindPFlag("verbosity", pflags.Lookup("verbosity"))

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// infoCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// infoCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
}
