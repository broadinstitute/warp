package error

import "fmt"

func HandleError(err error) {
	fmt.Printf("ERROR - %s", err.Error())
}
