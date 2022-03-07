package cmd

import (
	"fmt"
	"os"

	e "github.com/broadinstitute/warp/wreleaser/pkg/error"
	"github.com/broadinstitute/warp/wreleaser/pkg/releases"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

var version string

// listCmd represents the list command
var listCmd = &cobra.Command{
	Use:   "list",
	Short: "Display release information for a list of pipelines",
	Long: `'list' command expects one or more pipelines as arguments in the form:

        'wreleaser info list pipeline1 pipeline2 pipeline3 [flags]...'

        ** Note: If version flag is specified 'list' command will only accept one pipeline as argument **


Usage examples:

        'wreleaser info list Arrays' (display all releases for the Arrays pipeline)

        'wreleaser info list Arrays --version=v2.3.1' (display release information for Arrays pipeline v2.3.1)

        'wreleaser info list Arrays ExomeGermlineSingleSample' (display all releases for Arrays and ExomeGermlineSingleSample pipelines)

        'wreleaser info list CEMBA Optimus --latest' (display only the latest release for CEMBA and Optimus pipelines)`,

	Run: func(cmd *cobra.Command, args []string) {
		if len(args) == 0 {
			fmt.Fprintf(os.Stderr, "ERROR! Incompatible arguments - %v - 'list command expects at least one argument \n", args)
			fmt.Fprintf(os.Stderr, "Run 'wreleaser info list --help' to see example commands \n")
			os.Exit(1)
		}
		// Validate that args are valid pipelines
		releases.ValidateArgs(args)

		// Get raw list
		resp, err := releases.NewReleaseList()
		if err != nil {
			e.HandleError(err)
		}

		// Format the list and reduce to requested pipelines
		formatted, err := resp.FormatList(args)
		if err != nil {
			e.HandleError(err)
		}

		// Reduce the list to latest if requested
		if latest {
			// Verify user hasn't specified a version along with latest
			if version == "" {
				formatted.GetLatest()
				formatted.Print()
				os.Exit(0)
			} else {
				fmt.Fprintf(os.Stderr, "ERROR! Incompatible flags - 'latest'=%t   'version'=%s - Can't query unique VERSION and LATEST  - \n", latest, version)
				fmt.Fprintf(os.Stderr, "Run 'wreleaser info list --help' to see example commands \n")
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
				fmt.Fprintf(os.Stderr, "ERROR! Incompatible arguments - 'version'=%s    'pipeline args'=%v -  Can't query unique VERSION for multiple pipelines, please choose one pipeline \n", version, args)
				fmt.Fprintf(os.Stderr, "Run 'wreleaser info list --help' to see example commands \n")
				os.Exit(1)
			}
		}

		formatted.Print()
	},
}

func init() {
	infoCmd.AddCommand(listCmd)

	pflags := listCmd.PersistentFlags()

	pflags.StringVar(&version, "version", "", "Specific version number of pipeline to retrieve (can only use for single pipeline arguments)")

	viper.BindPFlag("version", pflags.Lookup("version"))
}
