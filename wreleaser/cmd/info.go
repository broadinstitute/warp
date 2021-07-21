package cmd

import (
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

var latest bool

// infoCmd represents the info command
var infoCmd = &cobra.Command{
	Use:   "info",
	Short: "Display release information for specified pipeline(s)",
	Long: `'info' command calls target subcommands

See subcommands:

        - all

        Example: 'wreleaser info all --help'

        - list

        Example 'wreleaser info list --help'

Default information (non-verbose) provided for the pipeline:
	- Release Name
	- Release Date
	- Release Notes
	- Release URL`,
}

func init() {
	rootCmd.AddCommand(infoCmd)

	pflags := infoCmd.PersistentFlags()

	pflags.BoolVarP(&latest, "latest", "l", false, "Retrieve only the latest releases for each pipeline")

	viper.BindPFlag("latest", pflags.Lookup("latest"))
}
