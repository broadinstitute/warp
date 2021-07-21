package error

import (
	"fmt"
	"os"
)

// HandleError prints errors
func HandleError(err error) {
	fmt.Fprintf(os.Stderr, "ERROR - %s \n", err.Error())
}
