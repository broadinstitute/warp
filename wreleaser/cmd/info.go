package cmd

import (
	"fmt"
	"log"
	"os"

	"github.com/spf13/cobra"
	"github.com/spf13/viper"

	"github.com/broadinstitute/warp/wreleaser/pkg/releases"
)

var (
	latest  bool
	version string
)

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
		releases.ValidateArgs(args)

		resp, err := releases.NewReleaseList()
		if err != nil {
			fmt.Print(err)
		}

		// Reduce the list to only requested pipelines
		formatted, err := resp.FormatList(args)
		if err != nil {
			fmt.Print(err)
		}

		// Reduce the list to latest if requested
		if latest {
			// Verify user hasn't specified a version along with latest
			if version == "" {
				formatted.GetLatest()
				formatted.Print()
				os.Exit(0)
			} else {
				log.Printf("ERROR! Incompatible arguments - 'latest'=%t   'version'=%s - Can't query unique VERSION and LATEST  -     \n", latest, version)
				log.Printf("Run 'wreleaser info --help' to see example commands")
				os.Exit(1)
			}
		}

		// Reduce the list to a specific version if requested
		if version != "" {
			// Verify user specified only one pipeline
			if len(args) == 1 {
				formatted.GetVersion(version, args[0])
				formatted.Print()
				os.Exit(0)
			} else {
				log.Printf("ERROR! Incompatible arguments - 'version'=%s    'pipeline args'=%v -  Can't query unique VERSION for multiple pipelines, please choose one pipeline    \n", version, args)
				log.Printf("Run 'wreleaser info --help' to see example commands")
				os.Exit(1)
			}
		}

		formatted.Print()
	},
}

func init() {
	rootCmd.AddCommand(infoCmd)

	pflags := infoCmd.PersistentFlags()

	pflags.BoolVarP(&latest, "latest", "l", false, "Retrieve only the latest releases for each pipeline")
	pflags.StringVar(&version, "version", "", "Version of the pipeline to query for")

	viper.BindPFlag("latest", pflags.Lookup("latest"))
	viper.BindPFlag("version", pflags.Lookup("version"))
}
