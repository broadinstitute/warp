package releases

import (
	"encoding/json"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"strings"

	"github.com/go-resty/resty/v2"

	e "github.com/broadinstitute/warp/wreleaser/pkg/error"
)

var (
	pageNumber   = 1
	responseSize = 100
	url          = "https://api.github.com/repos/broadinstitute/warp/releases"
)

// Release describes a Github release for a WARP pipeline
type Release struct {
	Name        string `json:"name"`
	ID          int    `json:"id"`
	PreRelease  bool   `json:"prerelease"`
	PublishedAt string `json:"published_at"`
	URL         string `json:"url"`
	HTMLURL     string `json:"html_url"`
	AssetsURL   string `json:"assets_url"`
	TarballURL  string `json:"tarball_url"`
	ZipballURL  string `json:"zipball_url"`
	Body        string `json:"body"`
}

// ReleaseList contains an array of Releases
type ReleaseList []Release

// NewReleaseList returns the full list of releases for all WARP pipelines
func NewReleaseList() (*ReleaseList, error) {
	// Check if the releases are cached
	// If not then we must fetch from Github
	if !cacheExists() {
		return makeNewList()
	}
	return listFromCache()
}

// FormatList formats a ReleaseList to a map[string]ReleaseList and removes prereleases
func (r *ReleaseList) FormatList() (*map[string]ReleaseList, error) {
	returnList := make(map[string]ReleaseList)

	for _, release := range *r {
		pipelineName := strings.Split(release.Name, "_")[0]

		if !release.PreRelease {
			if _, ok := returnList[pipelineName]; ok {
				returnList[pipelineName] = append(returnList[pipelineName], release)
			} else {
				returnList[pipelineName] = make(ReleaseList, 0)
				returnList[pipelineName] = append(returnList[pipelineName], release)
			}
		}
	}
	return &returnList, nil
}

// makeNewList creates the cache file and returns its values in a *ReleaseList
func makeNewList() (*ReleaseList, error) {
	var fullList ReleaseList
	cache, err := makeCache()
	if err != nil {
		e.HandleError(err)
	}
	defer cache.Close()

	client := resty.New()

	// Releases API chunks pages by 100
	// If responseSize is < 100 then we know we have hit the last page
	for responseSize == 100 {
		var temp ReleaseList

		resp, err := client.R().
			SetQueryParams(map[string]string{
				"per_page": fmt.Sprint(responseSize),
				"page":     fmt.Sprint(pageNumber),
			}).
			SetHeader("Accept", "application/vnd.github.v3+json").
			Get(url)
		if err != nil {
			return nil, err
		}

		if err := json.Unmarshal(resp.Body(), &temp); err != nil {
			return nil, err
		}

		fullList = append(fullList, temp...)

		responseSize = len(temp)

		pageNumber++
	}

	// Format and marshal full release list
	prettyJSON, err := json.MarshalIndent(fullList, "", "  ")
	if err != nil {
		return nil, err
	}

	// Write full formatted list to cache
	_, err = io.WriteString(cache, string(prettyJSON))
	if err != nil {
		return nil, err
	}

	// After creating the cachefile we can unmarshal from there and return
	return listFromCache()
}

// listFromCache returns a *ReleaseList based on the info in the cache file
func listFromCache() (*ReleaseList, error) {
	var (
		returnList ReleaseList
		cache      = getCacheDir()
	)

	cacheFile, err := os.Open(cache)
	if err != nil {
		return nil, err
	}
	defer cacheFile.Close()

	bytes, _ := ioutil.ReadAll(cacheFile)

	if err := json.Unmarshal(bytes, &returnList); err != nil {
		return nil, err
	}

	return &returnList, nil
}
