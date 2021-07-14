package cmd

import (
	"fmt"
	"os"

	e "github.com/broadinstitute/warp/wreleaser/pkg/error"
	"github.com/broadinstitute/warp/wreleaser/pkg/releases"
	"github.com/spf13/cobra"
)

// allCmd represents the all command
var allCmd = &cobra.Command{
	Use:   "all",
	Short: "Display release information for all pipelines",
	Long: `'all' command takes no arguments in the form:

        'wreleaser info all [flags]...'

Usage examples:

        'wreleaser info all' (display all releases for all pipelines)

        'wrleaser info all --latest (display latest releases for all pipelines)`,

	Run: func(cmd *cobra.Command, args []string) {
		if len(args) != 0 {
			fmt.Fprintf(os.Stderr, "ERROR! Incompatible arguments - %v - 'all' command takes no arguments \n", args)
			fmt.Fprintf(os.Stderr, "Run 'wreleaser info all -- help to see example commands \n")
			os.Exit(1)
		}

		// Get raw list
		resp, err := releases.NewReleaseList()
		if err != nil {
			e.HandleError(err)
		}

		// Format the list and reduce to requested pipelines (in this case all pipelines)
		formatted, err := resp.FormatList(args)
		if err != nil {
			e.HandleError(err)
		}

		// Reduce the list to latest release if requested
		if latest {
			formatted.GetLatest()
		}

		formatted.Print()
	},
}

func init() {
	infoCmd.AddCommand(allCmd)
}
