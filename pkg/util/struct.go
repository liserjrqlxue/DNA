package util

import (
	"errors"
	"fmt"
	"log"
	"log/slog"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/liserjrqlxue/goUtil/fmtUtil"
	mathUtil "github.com/liserjrqlxue/goUtil/math"
	"github.com/liserjrqlxue/goUtil/osUtil"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
)

type Seq struct {
	// Segment
	Start                   int
	End                     int
	TailOverlap             *Primer
	CohesiveTerminus        []byte
	SegmentCohesiveTerminus [][]byte
	SimplePlusMaxSegment    int

	Name    string
	RawSeq  string
	Seq     string
	Message string

	Length int

	Kmer   int
	KmerGC []float64
	// 1
	LowGC [][3]int
	// 2
	HighGC [][3]int

	// 3
	STR []*Feature
	// 4
	Poly []*Feature
	// 5
	Repeat []*Feature
	//  3+4+5
	TailRepeat []*Feature
	// 3+4
	MiddRepeat []*Feature

	SegmentPairs       []*PrimerPair
	CurrentSegmentPair *PrimerPair
	TiledPanel         [][8]*PrimerPair
	SplicerPanels      [][4][4]*PrimerPair

	// 临时变量
	start  int
	end    int
	offset int

	PrimerPairCount int // count of SegmentPairs

	UniversalUpperPrimer string
	UniversalLowerPrimer string

	RawSeqPrimer    *Primer
	PrimerS         *Primer
	PrimerSS        *Primer
	PrimerAS        *Primer
	PrimerASS       *Primer
	SegmentPrimers  []*Primer
	PrimerFs        []*Primer
	CandidatePrimer []*Primer
	CapturePrimers  []*Primer

	Note []string

	// output
	OutputPrefix string

	SegmentationTxt     *os.File
	FastaTxt            *os.File
	RepeatTxt           *os.File
	RepeatMiddTxt       *os.File
	RepeatMergedTxt     *os.File
	ExtraPrimerTxt      *os.File
	CandidateTxt        *os.File
	CandidateSegmentTxt *os.File
	SeqTxt              *os.File
	ThreePotTxt         *os.File
	NoteTxt             *os.File

	// segment split
	HomologousRecombination    bool
	SegmentPrimerRange         [2]int
	SegmentPrimerTmSuggest     float64
	SegmentPrimerLengthSuggest int
	// 候选分段overlap
	CandidateSegmentPrimer []*Primer
	Segments               []*Seq
	SegmentOverlapPrimer   *Primer
	SegmentLeftLimit       int
	SegmentRightLimit      int

	HardSwitch bool // hard模式启用开关
	Hard       bool
	Break      bool
}

func NewSeq(name, rawSeq, prefix string, hr, hardSwitch bool) *Seq {
	var s = &Seq{
		Name:         name,
		RawSeq:       rawSeq,
		OutputPrefix: prefix,

		RawSeqPrimer: NewPrimer(name, rawSeq, 0, len(rawSeq)),

		HomologousRecombination: hr,

		HardSwitch: hardSwitch,
	}
	s.Init()
	return s
}

func (s *Seq) Init() {
	if !ACGTN.MatchString(s.RawSeq) {
		slog.Error("Error (Invalid sequence)", "s.RawSeq", s.RawSeq)
		log.Fatalf("Error (Invalid sequence): [%s]", s.RawSeq)
	}

	s.Kmer = Kmer

	s.UniversalUpperPrimer = UniversalUpperPrimer
	s.UniversalLowerPrimer = UniversalLowerPrimer

	s.Seq = UniversalUpperPrimer + s.RawSeq + UniversalLowerPrimer
	s.Length = len(s.Seq)
	s.SegmentLeftLimit = len(UniversalUpperPrimer)
	s.SegmentRightLimit = s.Length - len(UniversalLowerPrimer)

	s.end = Min(s.start+SegPairLength, s.Length)
	s.offset = PrimerLengthMin

	if s.HomologousRecombination {
		s.SegmentPrimerRange = HomologousRecombinationSegmentPrimerRange
		s.SegmentPrimerLengthSuggest = 20
		s.SegmentPrimerTmSuggest = 58.0
	} else {
		s.SegmentPrimerRange = PrimerRange
		s.SegmentPrimerLengthSuggest = 30
		s.SegmentPrimerTmSuggest = TmRange[2]
	}
}

func (s *Seq) NewNoUni(name, sequence string) *Seq {
	if !ACGT.MatchString(sequence) {
		slog.Error("(Invalid sequence)", "sequence", sequence)
		log.Fatalf("Error (Invalid sequence): [%s]", sequence)
	}
	s.Name = name
	s.Seq = sequence
	s.Length = len(s.Seq)
	s.Kmer = Kmer

	s.end = Min(s.start+SegPairLength, s.Length)
	s.offset = PrimerLengthMin

	return s
}

func (s *Seq) CreateFiles() {
	s.FastaTxt = osUtil.Create(s.OutputPrefix + ".fa")
	s.RepeatTxt = osUtil.Create(s.OutputPrefix + ".repeat.bed")
	s.RepeatMiddTxt = osUtil.Create(s.OutputPrefix + ".repeat.midd.bed")
	s.RepeatMergedTxt = osUtil.Create(s.OutputPrefix + ".repeat.merged.bed")
	s.CandidateTxt = osUtil.Create(s.OutputPrefix + ".candidate.txt")
	s.CandidateSegmentTxt = osUtil.Create(s.OutputPrefix + ".candidate.Segment.txt")
	s.SeqTxt = osUtil.Create(s.OutputPrefix + ".seq")
}

func (s *Seq) CreateFiles4Step3() {
	s.SegmentationTxt = osUtil.Create(s.OutputPrefix + ".segmentation.txt")
	s.SeqTxt = osUtil.Create(s.OutputPrefix + ".seq")
	s.ThreePotTxt = osUtil.Create(s.OutputPrefix + ".3pot.txt")
	s.ExtraPrimerTxt = osUtil.Create(s.OutputPrefix + ".extra.txt")
	s.NoteTxt = osUtil.Create(s.OutputPrefix + ".note.txt")
}

func (s *Seq) Calculator() {
	s.CalculatorKmerGC()
	s.CalculatorPoly()
	s.CalculatorSTR()
	s.CalculatorRepeat()
}

func (s *Seq) CalculatorKmerGC() {
	var (
		GC     [SeqLengthMax]float64
		gcType []int
	)
	for i, c := range s.Seq {
		if c == 'C' || c == 'G' {
			for j := Max(0, i-s.Kmer); j <= i; j++ {
				GC[j]++
			}
		}
	}
	s.KmerGC = GC[:len(s.Seq)-s.Kmer+1]
	for i, f := range s.KmerGC {
		s.KmerGC[i] = f * 100 / float64(s.Kmer)
		if s.KmerGC[i] < gcRange[0] {
			gcType = append(gcType, -1)
		} else if s.KmerGC[i] > gcRange[1] {
			gcType = append(gcType, 1)
		} else {
			gcType = append(gcType, 0)
		}
	}
	var (
		lowReg  = [3]int{-1, -1, 1}
		highReg = [3]int{-1, -1, 2}
	)
	for i, t := range gcType {
		switch t {
		case 0:
			if lowReg[1] != -1 {
				s.LowGC = append(s.LowGC, lowReg)
				lowReg[0] = -1
				lowReg[1] = -1
			}
			if highReg[1] != -1 {
				s.HighGC = append(s.HighGC, highReg)
				highReg[0] = -1
				highReg[1] = -1
			}
		case -1:
			if lowReg[0] == -1 {
				lowReg[0] = i
			}
			lowReg[1] = i + s.Kmer

			if highReg[1] != -1 {
				s.HighGC = append(s.HighGC, highReg)
				highReg[0] = -1
				highReg[1] = -1
			}
		case 1:
			if highReg[0] == -1 {
				highReg[0] = i
			}
			highReg[1] = i + s.Kmer

			if lowReg[1] != -1 {
				s.LowGC = append(s.LowGC, lowReg)
				lowReg[0] = -1
				lowReg[1] = -1
			}
		}
	}
	if lowReg[1] != -1 {
		s.LowGC = append(s.LowGC, lowReg)
	}
	if highReg[1] != -1 {
		s.HighGC = append(s.HighGC, highReg)
	}
}

func (s *Seq) CalculatorSTR() {
	// 3 nt
	for kmer := 3; kmer < 5; kmer++ {
		var (
			End1 = s.Length - 2*kmer + 1
			End2 = s.Length - kmer + 1
		)
		for i := 0; i < End1; i++ {
			var (
				start = i
				end   = i
			)
			for j := end; j < End2; j += kmer {
				if s.Seq[j:j+kmer] != s.Seq[i:i+kmer] {
					break
				} else {
					end = j
				}
			}
			if end-start >= kmer {
				var str = &Feature{
					chr:   Name,
					start: start,
					end:   end + kmer,
					name:  fmt.Sprintf("STR:%s:%d", s.Seq[start:end+kmer], end-start+kmer),
				}
				s.STR = append(s.STR, str)
				if Verbose > 0 {
					slog.Info(str.String())
				}
			}
			i = end
		}
	}
}

// CalculatorPoly 至少 4 nt
func (s *Seq) CalculatorPoly() {
	for i := 0; i < s.Length-1; i++ {
		var (
			start = i
			end   = i
		)
		for j := i + 1; j < s.Length; j++ {
			if s.Seq[j] != s.Seq[i] {
				break
			} else {
				end = j
			}
		}
		if end-start >= 3 {
			var poly = &Feature{
				chr:   Name,
				start: start,
				end:   end + 1,
				name:  fmt.Sprintf("Poly:%s:%d", s.Seq[start:end+1], end-start+1),
			}
			s.Poly = append(s.Poly, poly)
			if Verbose > 0 {
				slog.Info("", "poly", poly)
			}
		}
		i = end
	}
}

func (s *Seq) CalculatorRepeat() {
	var (
		count      = make(map[string]int)
		kmer       = 5
		End        = s.Length - kmer + 1
		seq        = s.Seq
		repeatList []*Feature
		oldList    []*Feature
		allList    []*Feature
	)
	for i := 0; i < End; i++ {
		count[seq[i:i+kmer]]++
	}
	for i := 0; i < End; i++ {
		var str = seq[i : i+kmer]
		if count[str] > 1 {
			var repeat = &Feature{
				chr:   Name,
				start: i,
				end:   i + kmer,
				name:  fmt.Sprintf("Repeat:%s:%d", str, count[str]),
				seq:   str,
			}
			repeatList = append(repeatList, repeat)
		}
	}

	for kmer < 64 {
		kmer++
		End--

		if repeatList == nil {
			break
		}

		allList = append(allList, repeatList...)

		count = make(map[string]int)
		oldList = repeatList
		repeatList = []*Feature{}

		for _, reg := range oldList {
			if reg.end < s.Length {
				count[seq[reg.start:reg.end+1]]++
			}
		}
		for _, reg := range oldList {
			if reg.end < s.Length {
				var str = seq[reg.start : reg.end+1]
				if count[str] > 1 {
					var repeat = &Feature{
						chr:   Name,
						start: reg.start,
						end:   reg.end + 1,
						name:  fmt.Sprintf("Repeat:%s:%d", str, count[str]),
						seq:   str,
					}
					repeatList = append(repeatList, repeat)
				}
			}
		}
	}
	for _, repeat := range allList {
		if CheckGC70(repeat.seq) {
			s.Repeat = append(s.Repeat, repeat)
			if Verbose > 0 {
				slog.Info("", "repeat", repeat)
			}
		}
	}
	s.RepeatDeleteCovered()
}

func (s *Seq) RepeatDeleteCovered() {
	sort.Slice(s.Repeat, func(i, j int) bool {
		if s.Repeat[i].start == s.Repeat[j].start {
			return s.Repeat[i].end < s.Repeat[j].end
		}
		return s.Repeat[i].start < s.Repeat[j].start
	})
	var newRepeat []*Feature
	for i := 0; i < len(s.Repeat); i++ {
		var repeat1 = s.Repeat[i]
		var length = repeat1.end - repeat1.start
		var keep = true
		for j := 0; j < len(s.Repeat); j++ {
			var repeat2 = s.Repeat[j]
			if repeat2.start <= repeat1.start && repeat2.end >= repeat1.end && repeat2.end-repeat2.start > length {
				keep = false
				break
			}
		}
		if keep {
			newRepeat = append(newRepeat, repeat1)
		}
	}
	s.Repeat = newRepeat
}

func (s *Seq) CalTailRepeat() {
	var features []*Feature
	features = append(features, s.Poly...)
	features = append(features, s.STR...)
	features = append(features, s.Repeat...)
	sort.Slice(features, func(i, j int) bool {
		return features[i].start < features[j].start
	})

	features = MergeFeatures(s.Seq, features)
	s.TailRepeat = features
}

func (s *Seq) CalMiddRepeat() {
	var features []*Feature
	features = append(features, s.Poly...)
	features = append(features, s.STR...)
	sort.Slice(features, func(i, j int) bool {
		return features[i].start < features[j].start
	})
	features = MergeFeatures(s.Seq, features)
	s.MiddRepeat = features
}

func (s *Seq) FindSubSeq(start int) *Primer {
	var (
		primerCandidates []*Primer
		primerScore      []float64
	)
	for i := SubSeqRange[0]; i <= SubSeqRange[1]; i++ {
		if start+i > s.Length {
			continue
		}
		var primer = NewPrimer(s.Name+"-S", s.Seq, start, start+i)
		if primer.isInSubSeqRangeCandidate() {
			primerCandidates = append(primerCandidates, primer)
			primerScore = append(
				primerScore,
				math.Abs(SubSeqTmRange[2]-primer.Tm)*2+math.Abs(float64(primer.Length-45)),
			)
		}
	}
	if len(primerCandidates) > 0 {
		var _, index = FindMin(primerScore)
		return primerCandidates[index]
	}
	return nil
}

func (s *Seq) FindPrimers(start int) []*Primer {
	var (
		primerCandidates []*Primer
	)
	for i := PrimerRange[0]; i <= PrimerRange[1]; i++ {
		if start+i > s.Length {
			continue
		}
		var primer = NewPrimer(s.Name+"-S", s.Seq, start, start+i)
		if primer.Check(0, s.TailRepeat) {
			primerCandidates = append(primerCandidates, primer)
		}
	}
	return primerCandidates
}

// FindCapturePrimers
// ----------------------________  // start0 - end  : totalLength Tm <= 72
// __________------------________  // start  - end  : coreLength 20 delta 3 17-23, Tm = 56 extend = 6 56 - 62
// __________--------------------  // start  - end0 : totalLength Tm <=72
func (s *Seq) FindCapturePrimers(totalLength, rangeLength int, looseness bool) error {
	var (
		index = 1
		start = 0
		end   = start + rangeLength
	)
	for end <= s.Length {
		primer := s.Find1CapturePrimer(totalLength, start, end, index, looseness)
		if primer == nil {
			if end == s.Length {
				break
			}
			end = min(end+rangeLength, s.Length)
		} else {
			index++
			s.CapturePrimers = append(s.CapturePrimers, primer)
			start = end
			end = min(start+rangeLength, s.Length)
		}
	}
	if index == 1 {
		return errors.New("no capture primer found")
	}
	return nil
}

// coreLength = 20 delta = 3 // 17 - 23
// Tm = 56 extend = 6 56 - 62
// PairTm <= 72
// best in 200nt
// start0 -- start -- end -- end0
// [start, end)
// startL -- startR -- endL -- endR
// ------------------------________ // [startL , endL) : totalLength Tm <= 72
// __________---------------------- // [startR - endR) : totalLength Tm <=72
// __________--------------________ // [startR - endL) : coreLength 20 delta 3 17-23, Tm = 56 extend = 6 56 - 62
func (s *Seq) Find1CapturePrimer(totalLength, start, end, index int, looseness bool) *Primer {
	var (
		primerCandidates []*Primer
		primerScore      []float64
	)
	for startR := start; startR < end; startR++ {
		for coreLength := 17; coreLength <= 23; coreLength++ {
			endR := startR + totalLength
			endL := startR + coreLength
			startL := endL - totalLength

			if startL < 0 || endR > s.Length || endR > end {
				continue
			}

			var primer = NewPrimer(s.Name+"-C-"+strconv.Itoa(index), s.Seq, startR, startR+coreLength)
			if primer.Tm >= 56 && primer.Tm <= 62 && (looseness || primer.Check(0, s.TailRepeat)) {
				primerCandidates = append(primerCandidates, primer)
				primerScore = append(primerScore, math.Abs(primer.Tm-56)+math.Abs(float64(primer.Length-20))*0.6)
				slog.Info("Candidate Primer", "length", primer.Length, "Tm", primer.Tm)
			}
		}
	}

	candidatesCount := len(primerCandidates)
	slog.Debug("Find1CapturePrimer", "candidatesCount", candidatesCount)
	if len(primerCandidates) > 0 {
		// mathUtil.SortSlice2(primerCandidates, primerTmScore, primerLengthScore)
		mathUtil.SortSlice(primerCandidates, primerScore)
		for i, primer := range primerCandidates {
			left := NewPrimer(primer.Name+"-5F", s.Seq, primer.End-totalLength, primer.End)
			right := NewPrimer(primer.Name+"-3R", s.Seq, primer.Start, primer.Start+totalLength)
			// if left.Tm<=72 && right.Tm<=72 && left.Check(0, s.TailRepeat) && right.Check(0, s.TailRepeat) {
			slog.Debug("Candidate Primer", "i", i, "length", primer.Length, "extend", totalLength, "Tm", []float64{primer.Tm, left.Tm, right.Tm})
			if left.Tm <= 72 && right.Tm <= 72 {
				// 输出lef right Name Seq
				slog.Debug("find",
					"positon", primer.Start,
					// primer.Name, primer.Seq,
					left.Name, left.Seq,
					right.Name, right.Seq,
					"Tm", []float64{primer.Tm, left.Tm, right.Tm},
				)
				return primer
			}
		}
	}
	return nil
}

func (s *Seq) FindPrimerS() {
	var (
		primerCandidates []*Primer
		primerScore      []float64
	)
	for i := PrimerRange[0]; i <= PrimerRange[1]; i++ {
		if i > s.Length {
			continue
		}
		var primer = NewPrimer(s.Name+"-S", s.Seq, 20, i+20)
		if primer.Check(0, s.TailRepeat) {
			primerCandidates = append(primerCandidates, primer)
			primerScore = append(primerScore, math.Abs(TmRange[2]-primer.Tm))
		}
	}
	if len(primerCandidates) > 0 {
		var _, index = FindMin(primerScore)
		s.PrimerS = primerCandidates[index]
		s.PrimerSS = NewSimplePrimer(s.Name+"-SS", s.UniversalUpperPrimer+s.PrimerS.Seq)
	}
}

func (s *Seq) FindPrimerAS() {
	var (
		primerCandidates []*Primer
		primerScore      []float64
	)
	for i := PrimerRange[0]; i <= PrimerRange[1]; i++ {
		if i > s.Length-20 {
			continue
		}
		var primer = NewSimplePrimer(s.Name+"-AS", ReverseComplement(s.Seq[s.Length-i-20:s.Length-20]))
		if primer.Check(0, s.TailRepeat) {
			primerCandidates = append(primerCandidates, primer)
			primerScore = append(primerScore, math.Abs(TmRange[2]-primer.Tm))
		}
	}
	if len(primerCandidates) > 0 {
		var _, index = FindMin(primerScore)
		s.PrimerAS = primerCandidates[index]
		s.PrimerASS = NewSimplePrimer(s.Name+"-ASS", ReverseComplement(s.UniversalLowerPrimer)+s.PrimerAS.Seq)
	}
}

func (s *Seq) FindPrimer650(offset, count int) {
	var (
		primerCandidates []*Primer
		primerScore      []float64
	)
	// 使用尾引物
	if offset > s.Length-EdgeLength {
		return
	}
	for i := 0; i <= 50; i++ {
		var start = offset - i
		if start > s.Length-PrimerRange[0] {
			continue
		}
		for j := PrimerRange[0]; j <= PrimerRange[1]; j++ {
			if start+j > s.Length {
				continue
			}
			var primer = NewPrimer(
				s.Name+"-"+strconv.Itoa(count)+"F",
				s.Seq, start, start+j,
			)
			if primer.Check(0, s.TailRepeat) {
				primerCandidates = append(primerCandidates, primer)
				primerScore = append(primerScore, math.Abs(TmRange[2]-primer.Tm))
			}
		}
		start = offset + i
		if start > s.Length-PrimerRange[0] {
			continue
		}
		for j := PrimerRange[0]; j <= PrimerRange[1]; j++ {
			if start+j > s.Length {
				continue
			}
			var primer = NewPrimer(
				s.Name+"-"+strconv.Itoa(count)+"F",
				s.Seq, start, start+j,
			)
			if primer.Check(0, s.TailRepeat) {
				primerCandidates = append(primerCandidates, primer)
				primerScore = append(primerScore, math.Abs(TmRange[2]-primer.Tm))
			}
		}
	}

	if len(primerCandidates) > 0 {
		var _, index = FindMin(primerScore)
		s.SegmentPrimers = append(s.SegmentPrimers, primerCandidates[index])
	} else {
		s.SegmentPrimers = append(s.SegmentPrimers, NewSimplePrimer(s.Name+"-"+strconv.Itoa(count)+"F", ""))
	}
	s.FindPrimer650(offset+SangerLength, count+1)
}

var MaxLength = 1500

func (s *Seq) FindCandidateSubSeq() []*Primer {
	var primerCandidates []*Primer

	for start := 0; start < s.Length; start++ {
		for i := SubSeqRange[0]; i <= SubSeqRange[1]; i++ {
			if start+i > s.Length {
				continue
			}
			var p = NewPrimer(s.Name+"-S", s.Seq, start, start+i)
			if p.isInSubSeqRangeCandidate() {
				primerCandidates = append(primerCandidates, p)
			}
		}
	}
	return primerCandidates
}

// FindCandidateSegmentPrimer 寻找候选分节段primers
// FindCandidateSegmentPrimer finds candidate primers for a segment.
// It iterates over the sequence, starting from the beginning, and
// checks if each subsequence is a candidate primer. It appends the
// valid primers to the CandidateSegmentPrimer slice.
//
// The function checks if the subsequence length is within the range
// specified by SegmentPrimerRange. It also checks if the GC content is
// within the range (40, 68) if HomologousRecombination is true. It
// checks if the primer passes the GC zero condition and if it does
// not contain any local repeats within the TailRepeat.
func (s *Seq) FindCandidateSegmentPrimer() {
	var primerCandidates []*Primer

	// Iterate over the sequence
	for start := 0; start < s.Length; start++ {
		// Iterate over the subsequence lengths
		for i := s.SegmentPrimerRange[0]; i <= s.SegmentPrimerRange[1]; i++ {
			// Check if the subsequence is beyond the sequence length
			if start+i > s.Length {
				continue
			}
			// Create a new primer
			var p = NewPrimer(s.Name+"-S", s.Seq, start, start+i)
			// Check if the primer is a candidate
			if p.isInRangeCandidate() && p.isGCZero() && p.CheckLocalRepeat(0, s.TailRepeat) {
				// Check if homologous recombination is true
				if s.HomologousRecombination {
					// Check if the GC content is within the range
					if p.GC >= 40 && p.GC <= 68 {
						primerCandidates = append(primerCandidates, p)
					}
				} else {
					// Append the primer as a candidate
					primerCandidates = append(primerCandidates, p)
				}
			}
		}
	}
	// Set the candidate primers in the sequence
	s.CandidateSegmentPrimer = primerCandidates
}

func (s *Seq) FindCandidatePrimer() {
	var primerCandidates []*Primer

	for start := 0; start < s.Length; start++ {
		for i := PrimerRange[0]; i <= PrimerRange[1]; i++ {
			if start+i > s.Length {
				continue
			}
			var p = NewPrimer(s.Name+"-S", s.Seq, start, start+i)
			if p.isInRangeCandidate() && p.isGCZero() && p.CheckLocalRepeat(0, s.TailRepeat) {
				primerCandidates = append(primerCandidates, p)
			}
		}

	}
	s.CandidatePrimer = primerCandidates
}

func (s *Seq) FindCohesiveCandidatePrimer() {
	var primerCandidates []*Primer

	for start := 0; start < s.Length; start++ {
		for i := CohesiveTerminusSegmentPrimerRange[0]; i <= CohesiveTerminusSegmentPrimerRange[1]; i++ {
			if start+i > s.Length {
				continue
			}
			var p = NewPrimer(s.Name+"-S", s.Seq, start, start+i)
			if p.GC >= 40 && p.GC <= 68 && p.Tm >= 50 && p.Tm <= 70 && p.isGCZero() && p.CheckLocalRepeat(0, s.MiddRepeat) {
				primerCandidates = append(primerCandidates, p)
			}
		}
	}
	s.CandidateSegmentPrimer = primerCandidates
}

func (s *Seq) FindNextSegmentPair() bool {
	var pair *PrimerPair
	defer func() {
		if pair.Midd != nil {
			pair.SetLeftRight(s.Name, s.PrimerPairCount)

			s.SegmentPairs = append(s.SegmentPairs, pair)
			s.CurrentSegmentPair = pair
			s.offset = pair.Tail.Length
			s.start = pair.End - s.offset
			s.end = Min(s.start+SegPairLength, s.Length)
			s.PrimerPairCount++
		}
	}()

	var pairName = fmt.Sprintf("%s-%d", s.Name, s.PrimerPairCount+1)
	pair = s.NewPrimerPair(pairName, s.start, s.end)
	if pair.FindMiddTail(s.TailRepeat, s.MiddRepeat, false) {
		slog.Debug("FindMiddTail Done", "pair", pair.Name, "Length", pair.Length, "Start", pair.Start, "End", pair.End)
		slog.Debug("pair detail", "pair", pair.Name, "Head", pair.Head, "Tail", pair.Tail)
		return true
	}

	if s.HardSwitch {
		// Try Hard Mod
		slog.Debug("Try Hard", "s.Name", s.Name, "s.PrimerPairCount", s.PrimerPairCount)
		s.Hard = true
		s.end = Min(s.start+2*SegLengthHard-PrimerLengthMin, s.Length)
		pair = s.NewPrimerPair(pairName, s.start, s.end)
		if pair.FindMiddTail(s.TailRepeat, s.MiddRepeat, true) {
			slog.Debug("FindMiddTail Done", "pair", pair.Name, "Length", pair.Length, "Start", pair.Start, "End", pair.End)
			return true
		} else {
			slog.Debug("Can not find appropriate Midd+Tail Primer, try Segment", "pair.Seq", pair.Seq)
			s.Message = fmt.Sprintf("无法找到合适的引物对：[%d:%s]", pair.Start, pair.Seq)
			return false
		}
	}
	return false
}

func (s *Seq) CalAll() {
	s.Calculator()
	s.CalTailRepeat()
	s.CalMiddRepeat()
}

// write

// WriteFasta 序列 转换为 `FASTA` 格式
func (s *Seq) WriteFasta() {
	var fa = s.FastaTxt
	defer simpleUtil.DeferClose(fa)
	fmtUtil.Fprintf(fa, ">%s\n%s\n", s.Name, s.Seq)
}

func (s *Seq) WriteRepeat() {
	var out = s.RepeatTxt
	defer simpleUtil.DeferClose(out)

	for _, feature := range s.Poly {
		fmtUtil.Fprintf(out, "%s\n", feature.String())
	}
	for _, feature := range s.STR {
		fmtUtil.Fprintf(out, "%s\n", feature.String())
	}
	for _, feature := range s.Repeat {
		fmtUtil.Fprintf(out, "%s\n", feature.String())
	}
}

func (s *Seq) WriteMiddRepeat() {
	var out = s.RepeatMiddTxt
	defer simpleUtil.DeferClose(out)

	for _, feature := range s.Poly {
		fmtUtil.Fprintf(out, "%s\n", feature.String())
	}
	for _, feature := range s.STR {
		fmtUtil.Fprintf(out, "%s\n", feature.String())
	}
}

func (s *Seq) WriteTailRepeat() {
	var out = s.RepeatMergedTxt
	defer simpleUtil.DeferClose(out)

	for _, feature := range s.TailRepeat {
		fmtUtil.Fprintf(out, "%s\n", feature.String())
	}
}

// WritePrimerPairs 输出设计引物段结果，可作为 `SnapGene` `Import Primers from a list` 的输入
func (s *Seq) FindPrimerPairs() bool {
	for s.FindNextSegmentPair() {
		if s.CurrentSegmentPair.End == s.Length {
			return true
		}
	}
	// slog.Info("FindPrimerPairs fail, suggest fragementation", "s.Name", s.Name, "s.PrimerPairCount", s.PrimerPairCount)
	return false
}

func (s *Seq) WritePrimerPairs() {
	defer simpleUtil.DeferClose(s.SegmentationTxt)

	if len(s.Segments) > 0 {
		s.SegmentPairs = nil
		for _, segment := range s.Segments {
			s.SegmentPairs = append(s.SegmentPairs, segment.SegmentPairs...)
			// slog.Info("append(s.SplicerPanels, segment.SplicerPanels...)", "s.Name", s.Name, "segment.Name", segment.Name) //"segment.SplicerPanels", segment.SplicerPanels,

			s.SplicerPanels = append(s.SplicerPanels, segment.SplicerPanels...)
			s.Tile2Panel(segment.SegmentPairs)
		}
	} else {
		s.Segments = append(s.Segments, s)
		s.Tile2Panel(s.SegmentPairs)
	}
	for _, pair := range s.SegmentPairs {
		pair.SeqRef = s
		// log.Printf("%+v", pair)
		pair.Write()
	}
}

// 均分铺入8*2*6=96孔管
// 一个PrimerPair占用一对管孔
// 每个Segment独占两列管孔
func (s *Seq) Tile2Panel(SegmentPairs []*PrimerPair) {
	var (
		columnPair = int(math.Ceil(
			mathUtil.DivisionInt(len(SegmentPairs), 8),
		))
		columnCount = len(SegmentPairs) / columnPair
		columnExtra = len(SegmentPairs) % columnPair
		pairIndex   = 0
	)
	for i := 0; i < columnPair; i++ {
		var tile [8]*PrimerPair
		for j := 0; j < columnCount; j++ {
			tile[j] = SegmentPairs[pairIndex]
			pairIndex++
		}
		if columnExtra > 0 {
			tile[columnCount] = SegmentPairs[pairIndex]
			columnExtra--
			pairIndex++
		}
		s.TiledPanel = append(s.TiledPanel, tile)
	}
}

func (s *Seq) WriteSeq() {
	var out = s.SeqTxt
	defer simpleUtil.DeferClose(out)

	for _, primerPairs := range s.TiledPanel {
		for _, primerPair := range primerPairs {
			if primerPair != nil {
				if primerPair.Reverse {
					fmtUtil.Fprintf(out, "%s,%s\n", primerPair.Right.Name, ReverseComplement(primerPair.Right.Seq))
				} else {
					fmtUtil.Fprintf(out, "%s,%s\n", primerPair.Left.Name, primerPair.Left.Seq)
				}
			} else {
				fmtUtil.Fprintf(out, "%s,%s\n", "covering", "T")
			}
		}
		for _, primerPair := range primerPairs {
			if primerPair != nil {
				if primerPair.Reverse {
					fmtUtil.Fprintf(out, "%s,%s\n", primerPair.Left.Name, primerPair.Left.Seq)
				} else {
					fmtUtil.Fprintf(out, "%s,%s\n", primerPair.Right.Name, ReverseComplement(primerPair.Right.Seq))
				}
			} else {
				fmtUtil.Fprintf(out, "%s,%s\n", "covering", "T")
			}
		}
	}
}

func (s *Seq) Write3pot() {
	var out = s.ThreePotTxt
	defer simpleUtil.DeferClose(out)

	var (
		columnPair     = 3
		columnCount    = len(s.SegmentPairs) / columnPair
		columnExtra    = len(s.SegmentPairs) % columnPair
		pairIndex      = 0
		primerPairList []*PrimerPair
		seqList        [][2]string
	)

	for i := 0; i < columnPair; i++ {
		for j := 0; j < columnCount; j++ {
			pairIndex++
			primerPairList = append(primerPairList, s.SegmentPairs[pairIndex-1])
		}
		var index = columnCount
		if columnExtra > 0 {
			columnExtra--
			pairIndex++
			index++
			primerPairList = append(primerPairList, s.SegmentPairs[pairIndex-1])
		}
		for j := index; j < 8; j++ {
			primerPairList = append(primerPairList, nil)
		}
	}
	for i := 0; i < columnPair; i++ {
		for j := 0; j < 8; j++ {
			var index = i*8 + j
			var pair = primerPairList[index]
			if pair != nil {
				seqList = append(seqList, [2]string{
					primerPairList[index].Left.Name,
					primerPairList[index].Left.Seq,
				})
			} else {
				seqList = append(seqList, [2]string{
					"covering",
					"T",
				})
			}
		}
		for j := 0; j < 8; j++ {
			var index = i*8 + j
			var pair = primerPairList[index]
			if pair != nil {
				seqList = append(seqList, [2]string{
					primerPairList[index].Right.Name,
					ReverseComplement(primerPairList[index].Right.Seq),
				})
			} else {
				seqList = append(seqList, [2]string{
					"covering",
					"T",
				})
			}
		}
	}
	for _, seq := range seqList {
		fmtUtil.Fprintf(out, "%s,%s\n", seq[0], seq[1])
	}
}

func (s *Seq) PrintPrimer(out *os.File, primer *Primer) {
	if primer != nil {
		fmtUtil.Fprintf(out, "%s\t%s\t%s\t%d\t%.2f\n", primer.Name, primer.Seq, s.Name, primer.Length, primer.Tm)
	} else {
		fmtUtil.Fprintf(out, "%s-S\t%s\t%s\t\t\n", s.Name, "NULL", s.Name)
	}
}

func (s *Seq) PrintPrimerS(out *os.File) {
	var primer = s.PrimerS
	s.PrintPrimer(out, primer)
	if primer == nil {
		s.Note = append(s.Note, "No -S primer")
	}
}

func (s *Seq) PrintPrimerAS(out *os.File) {
	var primer = s.PrimerAS
	s.PrintPrimer(out, primer)
	if primer == nil {
		s.Note = append(s.Note, "No -AS primer")
	}
}

func (s *Seq) FindTailChangedPrimer() {
	s.FindPrimerS()
	s.FindPrimerAS()
	s.FindPrimer650(EdgeLength, 1)
}

func (s *Seq) WriteExtraPrimer() {
	var out = s.ExtraPrimerTxt
	defer simpleUtil.DeferClose(out)

	s.PrintPrimerS(out)

	if len(s.Segments) > 0 {
		for i, segment := range s.Segments {
			if i != 0 {
				segment.PrintPrimerS(out)
			}
			if i != len(s.Segments)-1 {
				segment.PrintPrimerAS(out)
			}
		}
	}

	s.PrintPrimerAS(out)

	// SegmentPrimers
	{
		for _, primer := range s.SegmentPrimers {
			s.PrintPrimer(out, primer)
		}
	}
}

func (s *Seq) WriteCandicatePrimer() {
	var out = s.CandidateTxt
	defer simpleUtil.DeferClose(s.CandidateTxt)

	s.FindCandidatePrimer()
	var candidates = s.CandidatePrimer

	for i, candidate := range candidates {
		fmtUtil.Fprintf(
			out,
			"%d\t%s\n",
			i, candidate.String("+"),
		)
	}
}

func (seq *Seq) WriteNote() {
	var out = seq.NoteTxt
	defer simpleUtil.DeferClose(out)

	var segPairCount = len(seq.SegmentPairs)
	if segPairCount > SegPairCountLimit {
		seq.Note = append(seq.Note, "Too many segment pairs")
		slog.Warn("Result primerPair", "segPairCount", segPairCount)
	} else {
		slog.Info("Result primerPair", "segPairCount", segPairCount)
	}

	if seq.Hard {
		slog.Warn("PrimerDesign result ", "Level", "Hard")
		seq.Note = append(seq.Note, "Hard Mode")
	} else {
		slog.Info("PrimerDesign result ", "Level", "Normal")
	}

	fmtUtil.Fprintf(out, "%s\t%s\n", seq.Name, strings.Join(seq.Note, ", "))
}

func (s *Seq) NewPrimerPair(name string, start, end int) *PrimerPair {
	var pair = &PrimerPair{
		SeqRef: s,
		Name:   name,
		Seq:    s.Seq[start:end],
		Start:  start,
		End:    end,
		Length: end - start,
	}
	pair.SetHead(s.offset)
	return pair
}

// 寻找分段序列
func (s *Seq) FindSegmentOverlapPrimer(n int) bool {
	var (
		start            = s.SegmentLeftLimit
		end              = s.SegmentRightLimit
		primerCandidates []*Primer
		primerScore      []float64
		pos              = (start*(n-1) + end) / n
		extend           = 30
		candidatePrimer  []*Primer
	)
	slog.Debug("FindSegmentOverlapPrimer", "len(s.CandidateSegmentPrimer)", len(s.CandidateSegmentPrimer))
	for _, primer := range s.CandidateSegmentPrimer {
		if primer.Start >= start && primer.End <= end {
			candidatePrimer = append(candidatePrimer, primer)
		}
	}
	if len(candidatePrimer) == 0 {
		slog.Debug("segment overlap primer not found", "start", start, "end", end)
		s.Break = true
		return false
	}
	s.SegmentOverlapPrimer = nil
	for s.SegmentOverlapPrimer == nil {
		for _, primer := range candidatePrimer {
			if primer.Start >= pos-extend && primer.End <= pos+extend {
				primerCandidates = append(primerCandidates, primer)
				primerScore = append(
					primerScore,
					math.Abs(primer.Tm-s.SegmentPrimerTmSuggest)*5+
						math.Abs(float64(primer.Length-s.SegmentPrimerLengthSuggest)),
				)
			}
		}
		if len(primerCandidates) > 0 {
			var _, index = FindMin(primerScore)
			s.SegmentOverlapPrimer = primerCandidates[index]
		} else {
			extend += 30
		}
	}
	slog.Debug("segment overlap primer got", "start", start, "end", end, "s.SegmentOverlapPrimer", s.SegmentOverlapPrimer.String("+"))
	return true
}

// 寻找分段序列
func (s *Seq) FindSegmentOverlapPrimer1(left, right, n int) bool {
	var (
		primerCandidates []*Primer
		primerScore      []float64
		pos              = left + (right-left)/n
		extend           = 30
		candidatePrimer  []*Primer
	)
	slog.Debug("", "len(s.CandidateSegmentPrimer)", len(s.CandidateSegmentPrimer))
	for _, primer := range s.CandidateSegmentPrimer {
		if primer.Start >= left && primer.End <= right {
			candidatePrimer = append(candidatePrimer, primer)
		}
	}
	if len(candidatePrimer) == 0 {
		slog.Debug("segment overlap primer not found", "left", left, "right", right)
		return false
	}
	s.SegmentOverlapPrimer = nil
	for s.SegmentOverlapPrimer == nil {
		for _, primer := range candidatePrimer {
			if primer.Start >= pos-extend && primer.End <= pos+extend {
				primerCandidates = append(primerCandidates, primer)
				primerScore = append(
					primerScore,
					math.Abs(primer.Tm-s.SegmentPrimerTmSuggest)*5+
						math.Abs(float64(primer.Length-s.SegmentPrimerLengthSuggest)),
				)
			}
		}
		if len(primerCandidates) > 0 {
			var _, index = FindMin(primerScore)
			s.SegmentOverlapPrimer = primerCandidates[index]
		} else {
			extend += 30
		}
	}
	slog.Debug("segment overlap primer got", "left", left, "right", right, "s.SegmentOverlapPrimer", s.SegmentOverlapPrimer.String("+"))

	return true
}

func (s *Seq) SegmentSplit(pics int) bool {
	s.FindCandidateSegmentPrimer()
	for _, p := range s.CandidateSegmentPrimer {
		fmtUtil.Fprintln(s.CandidateSegmentTxt, p.String("+"))
	}
	s.CandidateSegmentTxt.Close()

	if !s.SegmentSplit2() {
		slog.Warn("segment split 2 failed")
		s.Hard = false
		s.SegmentLeftLimit = len(UniversalUpperPrimer)
		s.SegmentRightLimit = s.Length - len(UniversalLowerPrimer)
		if !s.SegmentSplit3() {
			if s.Length < 2000 {
				slog.Error("segment split 3 failed, stop for too short sequence", "s.Name", s.Name, "s.Length", s.Length)
				return false
			}
			slog.Warn("segment split 3 failed")
			s.Hard = false
			s.SegmentLeftLimit = len(UniversalUpperPrimer)
			s.SegmentRightLimit = s.Length - len(UniversalLowerPrimer)
			if !s.SegmentSplit4() {
				slog.Warn("segment split 4 failed")
				s.Hard = false
				s.SegmentLeftLimit = len(UniversalUpperPrimer)
				s.SegmentRightLimit = s.Length - len(UniversalLowerPrimer)
				if !s.SegmentSplit5() {
					slog.Warn("segment split 5 failed")
					s.Hard = false
					s.SegmentLeftLimit = len(UniversalUpperPrimer)
					s.SegmentRightLimit = s.Length - len(UniversalLowerPrimer)
					if !s.SegmentSplit6() {
						slog.Warn("segment split 6 failed")
						s.Hard = false
						s.SegmentLeftLimit = len(UniversalUpperPrimer)
						s.SegmentRightLimit = s.Length - len(UniversalLowerPrimer)
						if !s.SegmentSplit7() {
							slog.Error("segment split 7 failed")
							return false
						}
						if !s.SegmentSplit7() {
							slog.Warn("segment split 7 failed")
							s.Hard = false
							s.SegmentLeftLimit = len(UniversalUpperPrimer)
							s.SegmentRightLimit = s.Length - len(UniversalLowerPrimer)
							if !s.SegmentSplitN(8) {
								slog.Error("segment split 8 failed")
								return false
							}
						}
					}
				}
			}
		}
	}
	for i, v := range s.Segments {
		s.Segments[i] = NewSeq(
			fmt.Sprintf("%s%c", s.Name, 'a'+i),
			v.RawSeq,
			fmt.Sprintf("%s%c", s.OutputPrefix, 'a'+i),
			s.HomologousRecombination,
			s.HardSwitch,
		)
		if !s.Segments[i].SimpleRun() {
			slog.Error("segment split failed", "i", i, "name", s.Segments[i].Name)
			return false
		}
	}
	return true
}

func (s *Seq) SegmentSplit2() bool {
	s.FindSegmentOverlapPrimer(2)
	slog.Debug(
		"Run segment split 2",
		"Name", s.Name,
		"Length", s.Length,
		"Range",
		fmt.Sprintf(
			"[%d-%d]:[%d-%d]",
			s.SegmentLeftLimit, s.SegmentRightLimit,
			s.SegmentOverlapPrimer.Start, s.SegmentOverlapPrimer.End,
		),
	)

	var seqa = NewSeq(
		s.Name+"a",
		s.Seq[len(UniversalUpperPrimer):s.SegmentOverlapPrimer.End],
		s.OutputPrefix+"a",
		s.HomologousRecombination,
		s.HardSwitch,
	)
	var oka = seqa.SimpleRun()
	var seqb = NewSeq(
		s.Name+"b",
		s.Seq[s.SegmentOverlapPrimer.Start:s.Length-len(UniversalLowerPrimer)],
		s.OutputPrefix+"b",
		s.HomologousRecombination,
		s.HardSwitch,
	)
	var okb = seqb.SimpleRun()
	if oka {
		slog.Info("Run segment split 1/2 success", "name", s.Name, "end", s.SegmentOverlapPrimer.End)
	} else {
		slog.Info("Run segment split 1/2 failed", "name", s.Name, "end", s.SegmentOverlapPrimer.End)
	}
	if okb {
		slog.Info("Run segment split 2/2 success", "name", s.Name, "start", s.SegmentOverlapPrimer.Start)
	} else {
		slog.Info("Run segment split 2/2 failed", "name", s.Name, "start", s.SegmentOverlapPrimer.Start)
	}
	if oka && okb {
		s.Hard = seqa.Hard && seqb.Hard
		s.Segments = []*Seq{seqa, seqb}
		return true
	} else {
		if oka {
			s.SegmentLeftLimit = Max(s.SegmentOverlapPrimer.Start, s.SegmentLeftLimit+1)
		}
		if okb {
			s.SegmentRightLimit = Min(s.SegmentOverlapPrimer.End, s.SegmentRightLimit-1)
		}
		if s.Break || !(oka || okb) {
			s.Break = false
			return false
		} else {
			if !oka && !okb {
				s.SegmentLeftLimit = Max(s.SegmentOverlapPrimer.Start, s.SegmentLeftLimit+1)
				s.SegmentRightLimit = Min(s.SegmentOverlapPrimer.End, s.SegmentRightLimit-1)
			}
			return s.SegmentSplit2()
		}
	}
}

func (s *Seq) SegmentSplit3() bool {
	s.FindSegmentOverlapPrimer(3)
	slog.Debug(
		"Run segment split 3",
		"Name", s.Name,
		"Length", s.Length,
		"Range",
		fmt.Sprintf(
			"[%d-%d]:[%d-%d]",
			s.SegmentLeftLimit, s.SegmentRightLimit,
			s.SegmentOverlapPrimer.Start, s.SegmentOverlapPrimer.End,
		),
	)

	var seqa = NewSeq(
		s.Name+"a",
		s.Seq[len(UniversalUpperPrimer):s.SegmentOverlapPrimer.End],
		s.OutputPrefix+"a",
		s.HomologousRecombination,
		s.HardSwitch,
	)
	var seqb = NewSeq(
		s.Name+"b",
		s.Seq[s.SegmentOverlapPrimer.Start:s.Length-len(UniversalLowerPrimer)],
		s.OutputPrefix+"b",
		s.HomologousRecombination,
		s.HardSwitch,
	)
	if seqa.SimpleRun() {
		slog.Info("Run segment split 1/3 success", "name", s.Name, "end", s.SegmentOverlapPrimer.End)
		slog.Debug(
			"Run segment split 2/3",
			"seqb.Name", seqb.Name,
			"seqa.Name", seqa.Name,
			"Range",
			fmt.Sprintf(
				"[%d-%d]",
				len(UniversalUpperPrimer),
				s.SegmentOverlapPrimer.End,
			),
		)
		seqb.FindCandidateSegmentPrimer()
		if seqb.SegmentSplit2() {
			s.Segments = append([]*Seq{seqa}, seqb.Segments...)
			return true
		} else if s.Break {
			s.Break = false
			return false
		} else {
			s.SegmentLeftLimit = Max(s.SegmentOverlapPrimer.Start, s.SegmentLeftLimit+1)
			return s.SegmentSplit3()
		}
	} else {
		slog.Info("Run segment split 1/3 failed", "name", s.Name, "end", s.SegmentOverlapPrimer.End)
		if s.Break {
			s.Break = false
			return false
		} else {
			s.SegmentRightLimit = Min(s.SegmentOverlapPrimer.End, s.SegmentRightLimit-1)
			return s.SegmentSplit3()
		}
	}
}

func (s *Seq) SegmentSplit4() bool {
	s.FindSegmentOverlapPrimer(4)
	slog.Debug(
		"Run segment split 4",
		"Name", s.Name,
		"Length", s.Length,
		"Range",
		fmt.Sprintf(
			"[%d-%d]:[%d-%d]",
			s.SegmentLeftLimit, s.SegmentRightLimit,
			s.SegmentOverlapPrimer.Start, s.SegmentOverlapPrimer.End,
		),
	)

	var seqa = NewSeq(
		s.Name+"a",
		s.Seq[len(UniversalUpperPrimer):s.SegmentOverlapPrimer.End],
		s.OutputPrefix+"a",
		s.HomologousRecombination,
		s.HardSwitch,
	)
	var seqb = NewSeq(
		s.Name+"b",
		s.Seq[s.SegmentOverlapPrimer.Start:s.Length-len(UniversalLowerPrimer)],
		s.OutputPrefix+"b",
		s.HomologousRecombination,
		s.HardSwitch,
	)
	if seqa.SimpleRun() {
		slog.Info("Run segment split 1/4 success", "name", s.Name, "end", s.SegmentOverlapPrimer.End)
		slog.Debug(
			"Run segment split 3/4",
			"seqb.Name", seqb.Name,
			"seqa.Name", seqa.Name,
			"Range",
			fmt.Sprintf(
				"[%d-%d]",
				len(UniversalUpperPrimer),
				s.SegmentOverlapPrimer.End,
			),
		)
		seqb.FindCandidateSegmentPrimer()
		if seqb.SegmentSplit3() {
			s.Segments = append([]*Seq{seqa}, seqb.Segments...)
			return true
		} else if s.Break {
			s.Break = false
			return false
		} else {
			s.SegmentLeftLimit = Max(s.SegmentOverlapPrimer.Start, s.SegmentLeftLimit+1)
			return s.SegmentSplit4()
		}
	} else {
		slog.Info("Run segment split 1/4 failed", "name", s.Name, "end", s.SegmentOverlapPrimer.End)
		if s.Break {
			s.Break = false
			return false
		} else {
			s.SegmentRightLimit = Min(s.SegmentOverlapPrimer.End, s.SegmentRightLimit-1)
			return s.SegmentSplit4()
		}
	}
}

func (s *Seq) SegmentSplit5() bool {
	s.FindSegmentOverlapPrimer(5)
	slog.Debug(
		"Run segment split 5",
		"Name", s.Name,
		"Length", s.Length,
		"Range",
		fmt.Sprintf(
			"[%d-%d]:[%d-%d]",
			s.SegmentLeftLimit, s.SegmentRightLimit,
			s.SegmentOverlapPrimer.Start, s.SegmentOverlapPrimer.End,
		),
	)

	var seqa = NewSeq(
		s.Name+"a",
		s.Seq[len(UniversalUpperPrimer):s.SegmentOverlapPrimer.End],
		s.OutputPrefix+"a",
		s.HomologousRecombination,
		s.HardSwitch,
	)
	var seqb = NewSeq(
		s.Name+"b",
		s.Seq[s.SegmentOverlapPrimer.Start:s.Length-len(UniversalLowerPrimer)],
		s.OutputPrefix+"b",
		s.HomologousRecombination,
		s.HardSwitch,
	)
	if seqa.SimpleRun() {
		slog.Info("Run segment split 1/5 success", "name", s.Name, "end", s.SegmentOverlapPrimer.End)
		slog.Debug(
			"Run segment split 4/5",
			"seqb.Name", seqb.Name,
			"seqa.Name", seqa.Name,
			"Range",
			fmt.Sprintf(
				"[%d-%d]",
				len(UniversalUpperPrimer),
				s.SegmentOverlapPrimer.End,
			),
		)
		seqb.FindCandidateSegmentPrimer()
		if seqb.SegmentSplit4() {
			s.Segments = append([]*Seq{seqa}, seqb.Segments...)
			return true
		} else if s.Break {
			s.Break = false
			return false
		} else {
			s.SegmentLeftLimit = Max(s.SegmentOverlapPrimer.Start, s.SegmentLeftLimit+1)
			return s.SegmentSplit5()
		}
	} else {
		slog.Info("Run segment split 1/5 failed", "name", s.Name, "end", s.SegmentOverlapPrimer.End)
		if s.Break {
			s.Break = false
			return false
		} else {
			s.SegmentRightLimit = Min(s.SegmentOverlapPrimer.End, s.SegmentRightLimit-1)
			return s.SegmentSplit5()
		}
	}
}

func (s *Seq) SegmentSplit6() bool {
	s.FindSegmentOverlapPrimer(6)
	slog.Debug(
		"Run segment split 6",
		"Name", s.Name,
		"Length", s.Length,
		"Range",
		fmt.Sprintf(
			"[%d-%d]:[%d-%d]",
			s.SegmentLeftLimit, s.SegmentRightLimit,
			s.SegmentOverlapPrimer.Start, s.SegmentOverlapPrimer.End,
		),
	)

	var seqa = NewSeq(
		s.Name+"a",
		s.Seq[len(UniversalUpperPrimer):s.SegmentOverlapPrimer.End],
		s.OutputPrefix+"a",
		s.HomologousRecombination,
		s.HardSwitch,
	)
	var seqb = NewSeq(
		s.Name+"b",
		s.Seq[s.SegmentOverlapPrimer.Start:s.Length-len(UniversalLowerPrimer)],
		s.OutputPrefix+"b",
		s.HomologousRecombination,
		s.HardSwitch,
	)
	if seqa.SimpleRun() {
		slog.Info("Run segment split 1/6 success", "name", s.Name, "end", s.SegmentOverlapPrimer.End)
		slog.Debug(
			"Run segment split 5/6",
			"seqb.Name", seqb.Name,
			"seqa.Name", seqa.Name,
			"Range",
			fmt.Sprintf(
				"[%d-%d]",
				len(UniversalUpperPrimer),
				s.SegmentOverlapPrimer.End,
			),
		)
		seqb.FindCandidateSegmentPrimer()
		if seqb.SegmentSplit5() {
			s.Segments = append([]*Seq{seqa}, seqb.Segments...)
			return true
		} else if s.Break {
			s.Break = false
			return false
		} else {
			s.SegmentLeftLimit = Max(s.SegmentOverlapPrimer.Start, s.SegmentLeftLimit+1)
			return s.SegmentSplit6()
		}
	} else {
		slog.Info("Run segment split 1/6 failed", "name", s.Name, "end", s.SegmentOverlapPrimer.End)
		if s.Break {
			s.Break = false
			return false
		} else {
			s.SegmentRightLimit = Min(s.SegmentOverlapPrimer.End, s.SegmentRightLimit-1)
			return s.SegmentSplit6()
		}
	}
}

func (s *Seq) SegmentSplit7() bool {
	s.FindSegmentOverlapPrimer(7)
	slog.Info(
		"Run segment split 7",
		"Name", s.Name,
		"Length", s.Length,
		"Range",
		fmt.Sprintf(
			"[%d-%d]:[%d-%d]",
			s.SegmentLeftLimit, s.SegmentRightLimit,
			s.SegmentOverlapPrimer.Start, s.SegmentOverlapPrimer.End,
		),
	)

	var seqa = NewSeq(
		s.Name+"a",
		s.Seq[len(UniversalUpperPrimer):s.SegmentOverlapPrimer.End],
		s.OutputPrefix+"a",
		s.HomologousRecombination,
		s.HardSwitch,
	)
	var seqb = NewSeq(
		s.Name+"b",
		s.Seq[s.SegmentOverlapPrimer.Start:s.Length-len(UniversalLowerPrimer)],
		s.OutputPrefix+"b",
		s.HomologousRecombination,
		s.HardSwitch,
	)
	if seqa.SimpleRun() {
		slog.Info("Run segment split 1/7 success", "name", s.Name, "end", s.SegmentOverlapPrimer.End)
		slog.Debug(
			"Run segment split 6/7",
			"seqb.Name", seqb.Name,
			"seqa.Name", seqa.Name,
			"Range",
			fmt.Sprintf(
				"[%d-%d]",
				len(UniversalUpperPrimer),
				s.SegmentOverlapPrimer.End,
			),
		)
		seqb.FindCandidateSegmentPrimer()
		if seqb.SegmentSplit6() {
			s.Segments = append([]*Seq{seqa}, seqb.Segments...)
			return true
		} else if s.Break {
			s.Break = false
			return false
		} else {
			s.SegmentLeftLimit = Max(s.SegmentOverlapPrimer.Start, s.SegmentLeftLimit+1)
			return s.SegmentSplit7()
		}
	} else {
		slog.Info("Run segment split 1/7 failed", "name", s.Name, "end", s.SegmentOverlapPrimer.End)
		if s.Break {
			s.Break = false
			return false
		} else {
			s.SegmentRightLimit = Min(s.SegmentOverlapPrimer.End, s.SegmentRightLimit-1)
			return s.SegmentSplit7()
		}
	}
}

func (s *Seq) SegmentSplitN(splitN int) bool {
	s.FindSegmentOverlapPrimer(splitN)
	slog.Info(
		"Run segment split", "splitN", splitN,
		"Name", s.Name,
		"Length", s.Length,
		"Range",
		fmt.Sprintf(
			"[%d-%d]:[%d-%d]",
			s.SegmentLeftLimit, s.SegmentRightLimit,
			s.SegmentOverlapPrimer.Start, s.SegmentOverlapPrimer.End,
		),
	)
	if splitN == 2 {
		return s.SegmentSplit2()
	}

	var seqa = NewSeq(
		s.Name+"a",
		s.Seq[len(UniversalUpperPrimer):s.SegmentOverlapPrimer.End],
		s.OutputPrefix+"a",
		s.HomologousRecombination,
		s.HardSwitch,
	)
	var seqb = NewSeq(
		s.Name+"b",
		s.Seq[s.SegmentOverlapPrimer.Start:s.Length-len(UniversalLowerPrimer)],
		s.OutputPrefix+"b",
		s.HomologousRecombination,
		s.HardSwitch,
	)
	if seqa.SimpleRun() {
		slog.Info("Run segment split 1 success", "splitN", splitN, "name", s.Name, "end", s.SegmentOverlapPrimer.End)
		slog.Debug(
			"Run segment split 2-N",
			"splitN", splitN,
			"seqb.Name", seqb.Name,
			"seqa.Name", seqa.Name,
			"Range",
			fmt.Sprintf(
				"[%d-%d]",
				len(UniversalUpperPrimer),
				s.SegmentOverlapPrimer.End,
			),
		)
		seqb.FindCandidateSegmentPrimer()
		if seqb.SegmentSplitN(splitN - 1) {
			s.Segments = append([]*Seq{seqa}, seqb.Segments...)
			return true
		} else if s.Break {
			s.Break = false
			return false
		} else {
			s.SegmentLeftLimit = Max(s.SegmentOverlapPrimer.Start, s.SegmentLeftLimit+1)
			return s.SegmentSplitN(splitN)
		}
	} else {
		slog.Info("Run segment split 1 failed", "splitN", splitN, "name", s.Name, "end", s.SegmentOverlapPrimer.End)
		if s.Break {
			s.Break = false
			return false
		} else {
			s.SegmentRightLimit = Min(s.SegmentOverlapPrimer.End, s.SegmentRightLimit-1)
			return s.SegmentSplitN(splitN)
		}
	}
}

func (s *Seq) recursiveSegmentSplit() {
	var maxSegmentNum = int(math.Ceil(float64(s.Length)/1000 + 2))
	s.Segments = make([]*Seq, maxSegmentNum)
	if !s.SegmentSplitNN(0, 2, maxSegmentNum, len(UniversalUpperPrimer), len(UniversalUpperPrimer), s.Length-len(UniversalLowerPrimer), s.Length-len(UniversalLowerPrimer)) {
		slog.Error("segment split failed", "s.Name", s.Name, "maxSegmentNum", maxSegmentNum)
		log.Fatalf("segment split %d failed:%s", maxSegmentNum, s.Name)
	}
}

func (s *Seq) SegmentSplitNN(currentN, totalN, maxN, start, left, right, end int) bool {
	// 范围 [left, right] 内找s.SegmentOverlapPrimer
	if !s.FindSegmentOverlapPrimer1(left, right, totalN-currentN) {
		if currentN == 0 && totalN < maxN {
			// 分更多段，起始
			return s.SegmentSplitNN(0, totalN+1, maxN, len(UniversalUpperPrimer), len(UniversalUpperPrimer), s.Length-len(UniversalLowerPrimer), s.Length-len(UniversalLowerPrimer))
		} else {
			return false
		}
	}
	var tag = string('a' + currentN)
	var seqCurrent = NewSeq(
		s.Name+tag,
		s.Seq[start:s.SegmentOverlapPrimer.End],
		s.OutputPrefix+tag,
		s.HomologousRecombination,
		s.HardSwitch,
	)
	if !seqCurrent.SimpleRun() {
		// 当前不合格，right左移
		right = Min(s.SegmentOverlapPrimer.End, right-1)
		return s.SegmentSplitNN(currentN, totalN, maxN, start, left, right, end)
	}

	// 递归处理下一段
	return s.SegmentSplitNN(currentN+1, totalN, maxN, start, left, right, end)
}

// 基础处理
func (seq *Seq) Step1() {
	slog.Debug("Run Step1", "seq.Name", seq.Name)
	seq.CalAll()
	seq.CreateFiles()
	seq.WriteFasta()
	seq.WriteRepeat()
	seq.WriteMiddRepeat()
	seq.WriteTailRepeat()
	seq.WriteCandicatePrimer()

	seq.FindTailChangedPrimer()
}

func (seq *Seq) Step3() {
	// Create File Handles
	seq.CreateFiles4Step3()
	// output after design
	seq.WritePrimerPairs()
	seq.WriteSeq()
	seq.Write3pot()
	seq.WriteExtraPrimer()
	seq.WriteNote()
}

func (seq *Seq) SimpleRun() bool {
	//InitLog(seq)
	seq.Step1()
	if seq.FindPrimerPairs() {
		if seq.PrimerPairCount > 16 {
			slog.Debug("Primer Pair Count > 16", "seq.Name", seq.Name, "seq.PrimerPairCount", seq.PrimerPairCount)
			return false
		} else {
			seq.BuildSplicerPanels()
		}
		// seq.Step3()
		return true
	}
	slog.Debug("SimpleRun Fail", "seq.Name", seq.Name)
	return false
}

func (seq *Seq) BuildSplicerPanels() {
	if seq.PrimerPairCount > 16 {
		slog.Error("Primer Pair Count > 16, Never BuildSplicerPanels", "seq.Name", seq.Name, "seq.PrimerPairCount", seq.PrimerPairCount)
		log.Fatal("Primer Pair Count > 16, Never BuildSplicerPanels")
		return
	}
	var SplicerPanel [4][4]*PrimerPair
	bk := SplitSegmentPair(seq.PrimerPairCount)
	for j := 0; j < 4; j++ {
		for i, pair := range seq.SegmentPairs[bk[j]:bk[j+1]] {
			SplicerPanel[j][i] = pair
			pair.Left.Name = fmt.Sprintf("%s_%d", seq.Name, 1+i+bk[j]+bk[j])
			pair.Right.Name = fmt.Sprintf("%s_%d", seq.Name, 1+i+bk[j]+bk[j+1])
		}
	}
	seq.SplicerPanels = [][4][4]*PrimerPair{SplicerPanel}
}

func Split2(a int, ceil bool) int {
	if ceil && a%2 == 1 {
		return a/2 + 1
	} else {
		return a / 2
	}
}

// 引物对拆分4份
func SplitSegmentPair(n int) [5]int {
	// N -> N1+N2 N1>=N2
	// 0..n -> 0..b-1 b--n -> 0:b b:n -> N[:b] N[b:]
	// N1: N[:b]
	// N2: N[b:]
	var b = Split2(n, true)
	// N1 -> N1a+N1b N1a>=N1b
	// 0..b-1 -> 0..b1-1 b--b -> 0:b1 b1:b -> N1[:b1] N1[b1:] -> N[:b1] N[b1:b]
	// N1a: N1[:b1] --> N[:b1]
	// N1b: N1[b1:] --> N[b1:b]
	var b1 = Split2(b, true)
	// N2 -> N2a+N1b N2a<=N2b
	// b..n -> b+0..b+b2-1 b+b2--n -> b:b+b2 b+b2:n -> N2[:b2] N1[b2:] -> N[b:b+b2] N[b+b2:n]
	// N1a: N2[:b2] --> N[b:b+b2]
	// N1b: N2[b2:] --> N[b+b2:]
	var b2 = Split2(n-b, false)
	return [5]int{0, b1, b, b + b2, n}
}

func SplitSegmentPairTest() {
	var data [16]int
	for i := 1; i <= 16; i++ {
		data[i-1] = i
		bk := SplitSegmentPair(i)
		fmt.Printf("%d:\t%+v\t%+v\t\n\tN1a:%+v\tN1b:%+v\tN2a:%+v\tN2b:%+v\n", i, data[:i], bk, data[:bk[0]], data[bk[0]:bk[1]], data[bk[1]:bk[2]], data[bk[2]:i])
	}
}
