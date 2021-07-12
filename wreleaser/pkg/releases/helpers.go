package releases

import (
	"os"
	"path/filepath"
)

func cacheExists(filename string) bool {
	info, err := os.Stat(filename)
	if os.IsNotExist(err) {
		return false
	}
	return !info.IsDir()
}

func makeCache(filename string) (*os.File, error) {
	if err := os.MkdirAll(filepath.Dir(filename), 0770); err != nil {
		return nil, err
	}
	return os.Create(filename)
}
