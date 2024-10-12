package main

/*
序列后面添加GC和Tm
*/

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log/slog"
	"math"
	"os"
	"runtime/pprof"
	"sort"
	"time"

	"github.com/liserjrqlxue/goUtil/osUtil"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
)

// cal GC content of a sequence

// flag
var (
	input = flag.String(
		"i",
		"",
		"input, one seq per line",
	)
	output = flag.String(
		"o",
		"",
		"output, one line output, format:\nSeq\tGCcontent\tTm",
	)
	cpuProfile = flag.String(
		"cpu",
		"",
		"write cpu profile to file",
	)
)

func main() {
	t0 := time.Now()
	flag.Parse()
	if *input == "" || *output == "" {
		flag.PrintDefaults()
		return
	}
	if *cpuProfile != "" {
		var LogCPUProfile = osUtil.Create(*cpuProfile)
		defer simpleUtil.DeferClose(LogCPUProfile)
		pprof.StartCPUProfile(LogCPUProfile)
		defer pprof.StopCPUProfile()
	}

	var (
		inF   *os.File
		outF  *os.File
		outGC *os.File
		outTm *os.File
		w     *bufio.Writer
		wGC   *bufio.Writer
		wTm   *bufio.Writer
	)

	// open input
	if *input == "-" {
		inF = os.Stdin
	} else {
		inF = osUtil.Open(*input)
		defer simpleUtil.DeferClose(inF)
	}

	// open output
	outF = osUtil.Create(*output)
	defer simpleUtil.DeferClose(outF)
	w = bufio.NewWriterSize(outF, 10*1024*1024)

	outGC = osUtil.Create(*output + ".gc")
	defer simpleUtil.DeferClose(outGC)
	wGC = bufio.NewWriter(outGC)

	outTm = osUtil.Create(*output + ".tm")
	defer simpleUtil.DeferClose(outTm)
	wTm = bufio.NewWriter(outTm)

	// 计数 GC
	var gcHist, tmHist = AddGC(inF, w)
	simpleUtil.CheckErr(writeHist(gcHist, wGC))
	simpleUtil.CheckErr(writeHist(tmHist, wTm))

	slog.Info("Done", "elapsed", time.Since(t0))
}

// Bytes2GC return GC content of a sequence round to 2 decimal
func Bytes2GC(seq []byte) float64 {
	gc := 0
	for _, c := range seq {
		switch c {
		case 'G', 'C':
			gc++
		}
	}
	return float64(gc) / float64(len(seq))
}

func Rount4(x float64) float64 {
	return math.Round(x * 10000)
}
func Rount2(x float64) float64 {
	return math.Round(x * 100)
}

// AddGC read io.Reader and write *bufio.Writer, add gc content, and return hist
func AddGC(in io.Reader, w *bufio.Writer) (gcHist, tmHist map[float64]int) {
	gcHist = make(map[float64]int)
	tmHist = make(map[float64]int)
	scan := bufio.NewScanner(in)
	for scan.Scan() {
		seq := scan.Bytes()
		gc := Bytes2GC(seq)
		gcContent := Rount4(gc) / 100
		gcHist[gcContent]++
		tm := CalculateTm(len(seq), gc)
		tmRound := Rount2(tm) / 100
		tmHist[tmRound]++
		fmt.Fprintf(w, "%s\t%.2f\t%.2f\n", seq, gcContent, tm)
	}
	simpleUtil.CheckErr(scan.Err())
	simpleUtil.CheckErr(w.Flush())
	return
}

const (
	A = 69.3 // 64.9
	B = 41.0
	C = 650.0 // 41*16.4=672.4
	// D = 16.6
	// Na = 0.05
)

func CalculateTm(length int, gc float64) float64 {
	return A + B*gc - C/float64(length)
}

// map2hist
func writeHist(hist map[float64]int, w *bufio.Writer) (err error) {
	var sortKey = make([]float64, 0, len(hist))
	for k := range hist {
		sortKey = append(sortKey, k)
	}
	sort.Float64s(sortKey)
	for _, k := range sortKey {
		fmt.Fprintf(w, "%.2f\t%d\n", k, hist[k])
	}
	err = w.Flush()
	return
}
