package releases

import (
	"encoding/json"
	"fmt"
	"log"
	"net/http"

	"github.com/spf13/viper"
)

var (
	client       http.Client
	pageNumber   = 1
	responseSize = 100
	URL          = "https://api.github.com/repos/broadinstitute/warp/releases?per_page=100&page="
)

// Release describes a github release
type Release struct {
	Url        string        `json:"url"`
	Id         string        `json:"id"`
	TagName    string        `json:"tag_name"`
	AssetsUrl  string        `json:"assets_url"`
	HtmlUrl    string        `json:"html_url"`
	PreRelease bool          `json:"prerelease"`
	Assets     []interface{} `json:"assets"`
	TarballUrl string        `json:"tarball_url"`
	ZipballUrl string        `json:"zipball_url"`
	Body       string        `json:"body"`
}

type ReleastList *[]Release

// NewReleaseList returns a new list of releases
func NewReleaseList() {

}

func Do() {
	//var allReleases ReleaseList
	cachedir := viper.GetString("cacheDir")
	// Check if the releases are cached, if not then fetch them
	if !cacheExists(cachedir) {
		cache, err := makeCache(cachedir)
		if err != nil {
			log.Println(err)
		}
		defer cache.Close()

		// Releases API sends 100 releases per page, loop through until responseSize is < 100
		for responseSize == 100 {
			tempList := &[]Release{}
			URL := URL + fmt.Sprint(pageNumber)

			req, err := http.NewRequest("GET", URL, nil)
			if err != nil {
				log.Println(err)
			}
			log.Printf(URL)

			req.Header = http.Header{
				"Accept": []string{"application/vnd.github.v3+json"},
			}

			res, err := client.Do(req)
			if err != nil {
				log.Println(err)
			}
			defer res.Body.Close()

			json.NewDecoder(res.Body).Decode(tempList)
			log.Println("hey")

			//allReleases.Releases = append(allReleases.Releases, temp.Releases...)
			responseSize = 10
			//pageNumber++

			//file, _ := json.MarshalIndent(temp.Releases, "", " ")
			fmt.Printf("%+v\n", tempList)
			//ioutil.WriteFile(cachedir, file)
		}

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
