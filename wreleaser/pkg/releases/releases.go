package releases

import (
	"encoding/json"
	"fmt"
	"net/http"
)

const (
	URL string = "https://api.github.com/repos/broadinstitute/warp/releases"
)

var client http.Client

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

func init() {

}

func Foo() {
	req, err := http.NewRequest("GET", URL, nil)
	if err != nil {
		fmt.Print(err)
	}

	req.Header = http.Header{
		"Accept": []string{"application/vnd.github.v3+json"},
	}

	res, err := client.Do(req)
	if err != nil {
		fmt.Print(err)
	}
	defer res.Body.Close()

	list := &[]Release{}

	json.NewDecoder(res.Body).Decode(list)

	fmt.Printf("%+v\n", (*list)[1])
}
