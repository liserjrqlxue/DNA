package main

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
		"input, two line for one Seq, format:\n\tName\n\tSeq",
	)
	output = flag.String(
		"o",
		"",
		"output, one line output, format:\n\tName GCcontent [Tm]\n\tSeq",
	)
	outputGC = flag.String(
		"gc",
		"",
		"output gc hist, format:\n\tGCcontent\tCount",
	)
	Tm = flag.Bool(
		"tm",
		false,
		"add tm",
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
		w     *bufio.Writer
		wGC   *bufio.Writer
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
	if *outputGC != "" {
		outGC = osUtil.Create(*outputGC)
		defer simpleUtil.DeferClose(outGC)
		wGC = bufio.NewWriter(outGC)
	}

	// 计数 GC
	var gcHist = simpleUtil.HandleError(AddGC(inF, w, *Tm))
	if wGC != nil {
		simpleUtil.CheckErr(writeHist(gcHist, wGC))
	}

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
	return math.Round(10000*float64(gc)/float64(len(seq))) / 100
}

// AddGC read io.Reader and write *bufio.Writer, add gc content, and return hist
func AddGC(in io.Reader, w *bufio.Writer, Tm bool) (hist map[float64]int, err error) {
	hist = make(map[float64]int)
	scan := bufio.NewScanner(in)
	for scan.Scan() {
		name := scan.Bytes()
		seq := scan.Bytes()
		gcContent := Bytes2GC(seq)
		hist[gcContent]++
		if Tm {
			tm := CalculateTm(len(seq), gcContent)
			fmt.Fprintf(w, "%s %.2f %.2f\n%s\n", name, gcContent, tm, seq)
		} else {
			fmt.Fprintf(w, "%s %.2f\n%s\n", name, gcContent, seq)
		}
	}
	err = scan.Err()
	if err != nil {
		return
	}
	err = w.Flush()
	if err != nil {
		return
	}
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
	return A + B*gc/100 - C/float64(length)
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
