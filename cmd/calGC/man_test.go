package main

import (
	"fmt"
	"strconv"
	"testing"
)

func BenchmarkSprintf(b *testing.B) {
	var f float64 = 123.456
	for i := 0; i < b.N; i++ {
		_ = fmt.Sprintf("%.2f", f)
	}
}

func BenchmarkFormatFloat(b *testing.B) {
	var f float64 = 123.456
	for i := 0; i < b.N; i++ {
		_ = strconv.FormatFloat(f, 'f', 2, 64)
	}
}
