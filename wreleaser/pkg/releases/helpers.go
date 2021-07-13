package releases

import (
	"os"
	"path/filepath"

	"github.com/spf13/viper"
)

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
