package error

import "fmt"

// HandleError prints errors
func HandleError(err error) {
	fmt.Printf("ERROR - %s", err.Error())
}
