package cmd

import (
	"encoding/json"
	"fmt"

	"github.com/spf13/cobra"

	"github.com/broadinstitute/warp/wreleaser/pkg/releases"
)

var latest bool

// infoCmd represents the info command
var infoCmd = &cobra.Command{
	Use:   "info",
	Short: "Display information for specified pipeline(s)",
	Long: `This command will display information for specified pipeline(s)

Default information (non-verbose) provided for the pipeline:
	- Release Name
	- Release Id
	- Release Date
	- Release Notes
	- Various URLs (Tarball, Zip, Assets, Raw Html)

Usage examples:

'wreleaser info' (display release information for all available pipelines)

'wreleaser info -l' (display only latest (most recent) release information for all available pipelines)

'wreleaser info CEMBA ExomeGermlineSample Optimus' (display info for the CEMBA, ExomeGermline and Optimus pipelines)

`,
	Run: func(cmd *cobra.Command, args []string) {
		resp, err := releases.NewReleaseList()
		if err != nil {
			fmt.Print(err)
		}

		formatted, err := resp.FormatList()
		if err != nil {
			fmt.Print(err)
		}

		// Format and marshal full release list
		prettyJSON, err := json.MarshalIndent(*formatted, "", "  ")
		if err != nil {
			fmt.Print(err)
		}
		fmt.Print(string(prettyJSON))
	},
}

func init() {
	rootCmd.AddCommand(infoCmd)

	pflags := infoCmd.PersistentFlags()

	pflags.BoolVarP(&latest, "latest", "l", false, "Retrieve only the latest releases for each pipeline")

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// infoCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// infoCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
}
