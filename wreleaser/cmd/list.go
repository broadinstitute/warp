package cmd

import (
	"fmt"
	"os"

	"github.com/spf13/cobra"
)

// listCmd represents the list command
var listCmd = &cobra.Command{
	Use:   "list",
	Short: "Display release information for a list of pipelines",
	Long: `'list' command expects one or more pipelines as arguments in the form:

        'wreleaser info list pipeline1 pipeline2 pipeline3 [flags]...'

Usage examples:

        'wreleaser info list Arrays' (display all releases for the Arrays pipeline)

        'wreleaser info list Arrays ExomeGermline' (display all releases for Arrays and ExomeGermline pipelines)

        'wreleaser info list CEMBA Optimus --latest' (display only the latest release for CEMBA and Optimus pipelines)`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Fprintln(os.Stderr)
	},
}

func init() {
	infoCmd.AddCommand(listCmd)

	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// listCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// listCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
}
