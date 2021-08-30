package cmd

import (
	"fmt"
	"os"

	"github.com/spf13/cobra"
	"github.com/spf13/viper"

	e "github.com/broadinstitute/warp/wreleaser/pkg/error"
)

// resetCmd represents the reset command
var resetCmd = &cobra.Command{
	Use:   "reset",
	Short: "Clear the directory which caches WARP releases",
	Long: `'reset' command clears the directory of cached WARP releases

Usage example:

        'wreleaser reset' (clear the default cache located at $HOME/.wreleaser/cache.json)

        'wreleaser reset --cachedir=/path/to/custom' (clear the custom cache located at /path/to/custom)`,

	Run: func(cmd *cobra.Command, args []string) {
		cache := viper.GetString("cachedir")

		if err := os.Remove(cache); err != nil {
			e.HandleError(err)
			os.Exit(1)
		}

		fmt.Printf("Cache %s deleted", cache)
	},
}

func init() {
	rootCmd.AddCommand(resetCmd)
}
