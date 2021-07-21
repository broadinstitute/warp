package releases

import (
	"encoding/json"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"regexp"
	"strings"

	"github.com/go-resty/resty/v2"
	"github.com/spf13/viper"

	e "github.com/broadinstitute/warp/wreleaser/pkg/error"
)

var (
	pageNumber   = 1
	responseSize = 100
	versionFound = false
	url          = "https://api.github.com/repos/broadinstitute/warp/releases"
)

// Release describes a Github release for a WARP pipeline
type Release struct {
	Name        string `json:"name"`
	PreRelease  *bool  `json:"prerelease,omitempty"`
	PublishedAt string `json:"published_at"`
	HTMLURL     string `json:"html_url"`
	Body        string `json:"body"`
}

// ReleaseList contains an array of Releases
type ReleaseList []Release

// ReleaseListFormatted maps each pipeline name to all of its releases
type ReleaseListFormatted map[string]ReleaseList

// NewReleaseList returns the full list of releases for all WARP pipelines
func NewReleaseList() (*ReleaseList, error) {
	// Check if the releases are cached
	// If not then we must fetch from Github
	if !cacheExists() {
		return makeNewList()
	}
	return listFromCache()
}

// FormatList formats a ReleaseList to a map[string]ReleaseList, removes prereleases and filters to requested pipelines
func (r *ReleaseList) FormatList(requestedPipelines []string) (*ReleaseListFormatted, error) {
	returnList := make(ReleaseListFormatted)

	for _, release := range *r {
		pipelineName := strings.Split(release.Name, "_")[0]

		if !*release.PreRelease {
			// remove PreRelease field from struct
			release.PreRelease = nil
			if _, ok := returnList[pipelineName]; ok {
				returnList[pipelineName] = append(returnList[pipelineName], release)
			} else {
				// If no arguments then add all pipelines, otherwise only add requested pipelines from args
				if len(requestedPipelines) == 0 || hasValue(requestedPipelines, pipelineName) {
					returnList[pipelineName] = make(ReleaseList, 0)
					returnList[pipelineName] = append(returnList[pipelineName], release)
				}
			}
		}
	}
	return &returnList, nil
}

// GetLatest modifies a ReleaseListFormatted to only include the latest releases
func (r *ReleaseListFormatted) GetLatest() {
	for key, element := range *r {
		(*r)[key] = []Release{element[0]}
	}
}

// GetVersion modifies a ReleaseListFormmated to only include the requested version
func (r *ReleaseListFormatted) GetVersion(version string, pipeline string) {
	for key, list := range *r {
		for _, release := range list {
			releaseVersion := strings.Split(release.Name, "_")[1]
			if version == releaseVersion {
				(*r)[key] = []Release{release}
				versionFound = true
			}
		}
	}
	if versionFound == false {
		fmt.Fprintf(os.Stderr, "ERROR - Unable to find' 'version'=%s for 'pipeline'=%s \n", version, pipeline)
		fmt.Fprintf(os.Stderr, "Run 'wreleaser info list %s' to see all release versions for the %s pipeline \n", pipeline, pipeline)
		os.Exit(1)
	}
}

// Print prints the ReleaseListFormmated to stdout
func (r *ReleaseListFormatted) Print() {
	prettyJSON, err := json.MarshalIndent(*r, "", "  ")
	if err != nil {
		e.HandleError(err)
	}

	output := viper.GetString("output")

	if output != "" {
		f, err := os.Create(output)
		if err != nil {
			e.HandleError(err)
		}

		// Write full formatted list to cache
		_, err = io.WriteString(f, string(prettyJSON))
		if err != nil {
			e.HandleError(err)
		}
	} else {
		fmt.Print(string(prettyJSON))
	}
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

		// replace all escape characters in response body with single space
		// replace all additional whitespace
		for i, release := range temp {
			tempReplace := strings.ReplaceAll(release.Body, "\r\n", " ")
			space := regexp.MustCompile(`\s+`)
			temp[i].Body = space.ReplaceAllString(tempReplace, " ")
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
