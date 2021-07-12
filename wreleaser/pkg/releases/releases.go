package releases

import (
	"encoding/json"
	"fmt"
	"log"
	"net/http"

	"github.com/go-resty/resty/v2"
	"github.com/spf13/viper"
)

var (
	client       http.Client
	pageNumber   = 1
	responseSize = "100"
	URL          = "https://api.github.com/repos/broadinstitute/warp/releases"
)

// Release describes a Github release for a WARP pipeline
type Release struct {
	Url        string        `json:"url"`
	Id         int           `json:"id"`
	TagName    string        `json:"tag_name"`
	AssetsUrl  string        `json:"assets_url"`
	HtmlUrl    string        `json:"html_url"`
	PreRelease bool          `json:"prerelease"`
	Assets     []interface{} `json:"assets"`
	TarballUrl string        `json:"tarball_url"`
	ZipballUrl string        `json:"zipball_url"`
	Body       string        `json:"body"`
}

type ReleaseList []Release

// NewReleaseList returns the full list of releases for all WARP pipelines
func NewReleaseList() {
	cachedir := viper.GetString("cacheDir")
	// Check if the releases are cached, if not then fetch them
	if !cacheExists(cachedir) {
		cache, err := makeCache(cachedir)
		if err != nil {
			log.Println(err)
		}
		defer cache.Close()

		var list ReleaseList
		client := resty.New()

		// Releases API sends 100 releases per page, loop through until responseSize is < 100
		resp, err := client.R().
			SetQueryParams(map[string]string{
				"per_page": responseSize,
				"page":     fmt.Sprint(pageNumber),
			}).
			SetHeader("Accept", "application/vnd.github.v3+json").
			Get(URL)
		if err != nil {
			log.Println(err)
		}
		log.Println(len(resp.Body()))
		if err := json.Unmarshal(resp.Body(), &list); err != nil {
			log.Println(err)
		}

		pretty, err := json.MarshalIndent(list[1], "", "  ")
		if err != nil {
			log.Println(err)
		}
		fmt.Printf("%s\n", pretty)

		//fmt.Printf("%+v", list[1])
		//log.Println(len(list))

	}

}

//func Foo() {
//	req, err := http.NewRequest("GET", URL, nil)
//	if err != nil {
//		log.Println(err)
//	}
//
//	req.Header = http.Header{
//		"Accept": []string{"application/vnd.github.v3+json"},
//	}
//
//	res, err := client.Do(req)
//	if err != nil {
//		fmt.Print(err)
//	}
//	defer res.Body.Close()
//
//	list := &[]Release{}
//
//	json.NewDecoder(res.Body).Decode(list)
//	fmt.Print(len((*list)))
//
//	//fmt.Printf("%+v\n", (*list)[1])
//}
