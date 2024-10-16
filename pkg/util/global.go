package util

import "regexp"

// const
const (
	seqLengthMax      = 1024 * 12
	Kmer              = 15
	PrimerLengthMin   = 20
	PrimerLength      = 20
	PrimerLengthMax   = 15
	SegmentExtend     = 10
	SegPairCountLimit = 16
	// 末端GC长度
	terminalGC = 2
	outRepeat  = 2
)

// limitation
var (
	SegLength     = 70
	SegLengthHard = 90
	SegPairLength = SegLength*2 - PrimerLengthMin

	// 引物长度范围
	PrimerRange = [2]int{20, 30}
	// 同源重组分段引物长度范围
	HomologousRecombinationSegmentPrimerRange = [2]int{15, 25}
	CohesiveTerminusSegmentPrimerRange        = [2]int{15, 20}

	SegLengthRange = [2]int{55, 90}
	TmRange        = []float64{58.0, 70.0, 62.0}
	gcRange        = []float64{40.0, 60.0, 50.0}
	SubSeqRange    = [2]int{41, 49}
	SubSeqTmRange  = []float64{72.0, 76.0, 74.0}
	SubSeqGcRange  = []float64{42.0, 50.0, 46.0}

	Name    string
	Verbose int
)

// regexp
var (
	// ACGT valid sequence
	ACGT  = regexp.MustCompile(`^[ACGT]*$`)
	ACGTN = regexp.MustCompile(`^[ACGTNK]*$`)
)

var (
	UniversalUpperPrimer string
	UniversalLowerPrimer string
)
