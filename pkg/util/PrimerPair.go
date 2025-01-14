package util

import (
	"fmt"
	"log/slog"
	"math"

	"github.com/liserjrqlxue/goUtil/fmtUtil"
)

type PrimerPair struct {
	Name string
	Seq  string

	Length int
	Start  int
	End    int

	Breakpoints [6]int

	Left    *Primer
	Right   *Primer
	Reverse bool

	Head *Primer
	Tail *Primer
	Midd *Primer

	SeqRef *Seq

	// CPs Chemical Primers
	CPs []*Primer
}

func (pair *PrimerPair) SetHead(offset int) {
	pair.Head = NewPrimer(
		pair.Name+"_1",
		pair.Seq,
		0,
		offset,
	)
}

func (pair *PrimerPair) SetTail(offset int) {
	pair.Tail = NewPrimer(
		pair.Name+"_3",
		pair.Seq,
		pair.Length-offset,
		pair.Length,
	)
}

// FindTail 终点固定，不定长
func (pair *PrimerPair) FindTail(repeat []*Feature, overwrite bool) {
	slog.Debug("FindTail", "pair", pair.Name, "Length", pair.Length, "Start", pair.Start, "End", pair.End)
	if pair.Tail != nil && !overwrite {
		return
	}
	var (
		primerCandidates []*Primer
		primerScore      []float64
	)
	for i := PrimerRange[0]; i <= PrimerRange[1]; i++ {
		if i > pair.Length {
			continue
		}
		var primer = NewPrimer(pair.Name+"_3", pair.Seq, pair.Length-i, pair.Length)
		if primer.Check(pair.Start, repeat) && primer.IsUnambiguous() {
			primerCandidates = append(primerCandidates, primer)
			primerScore = append(primerScore, math.Abs(TmRange[2]-primer.Tm))
		}
	}
	if len(primerCandidates) > 0 {
		var _, index = FindMin(primerScore)
		pair.Tail = primerCandidates[index]
	}
}

// FindMidd 范围内，不定长
func (pair *PrimerPair) FindMidd(repeat []*Feature, hard bool) {
	var (
		primerCandidates []*Primer
		primerScore      []float64
		segLength        = SegLength
	)

	for start := Max(0, pair.Length-segLength); start < segLength; start++ {
		for i := PrimerRange[0]; i <= PrimerRange[1]; i++ {
			if start+i > segLength || start+i > pair.Length {
				continue
			}
			var primer = NewPrimer(pair.Name+"_2", pair.Seq, start, start+i)
			if primer.Check(pair.Start, repeat) && primer.IsUnambiguous() {
				primerCandidates = append(primerCandidates, primer)
				primerScore = append(primerScore, math.Abs(TmRange[2]-primer.Tm))
			}
		}
	}
	if len(primerCandidates) == 0 && hard {
		for start := segLength; start < SegLengthHard; start++ {
			for i := PrimerRange[0]; i <= PrimerRange[1]; i++ {
				if start+i > segLength || start+i > pair.Length {
					continue
				}
				var primer = NewPrimer(pair.Name+"_2", pair.Seq, start, start+i)
				if primer.Check(pair.Start, repeat) && primer.IsUnambiguous() {
					primerCandidates = append(primerCandidates, primer)
					primerScore = append(primerScore, math.Abs(TmRange[2]-primer.Tm))
				}
			}
		}
	}
	if len(primerCandidates) > 0 {
		var _, index = FindMin(primerScore)
		pair.Midd = primerCandidates[index]
	}
}

func (pair *PrimerPair) RightCut() {
	pair.End--
	pair.Length--
	pair.Seq = pair.Seq[:pair.Length]
	pair.Tail = nil
}

// FindMiddTail 终点不固定，不定长
func (pair *PrimerPair) FindMiddTail(tailRepeat, middRepeat []*Feature, hard bool) bool {
	slog.Debug("FindMiddTail", "pair", pair.Name, "Length", pair.Length, "Start", pair.Start, "End", pair.End)
	// 到末端
	if pair.SeqRef.Length == pair.End {
		pair.SetTail(len(UniversalLowerPrimer))
		pair.FindMidd(pair.SeqRef.MiddRepeat, hard)
		if pair.Midd != nil {
			return true
		}
	}

	for pair.Length >= pair.Head.Length+SegmentExtend {
		pair.FindTail(tailRepeat, false)
		if pair.Tail != nil {
			pair.FindMidd(middRepeat, hard)
			if pair.Midd != nil {
				return true
			}
		}
		pair.RightCut()
	}
	if hard {
		slog.Debug("Can not find appropriate Midd+Tail Primer, try Segment", "pair.Seq", pair.Seq)
	} else {
		slog.Debug("Can not find appropriate Midd+Tail Primer, try Hard", "pair.Seq", pair.Seq)
	}
	return false
}

func (pair *PrimerPair) SetLeftRight(name string, index int) {
	pair.Left = NewPrimer(
		fmt.Sprintf("%s-%d", name, 2*index+1),
		pair.Seq, 0, pair.Midd.End,
	)
	pair.Right = NewPrimer(
		fmt.Sprintf("%s-%d", name, 2*index+2),
		pair.Seq, pair.Midd.Start, pair.Length,
	)
}

func (pair *PrimerPair) Write() {
	var (
		out   = pair.SeqRef.SegmentationTxt
		left  = pair.Left
		right = pair.Right
		Head  = pair.Head
		Midd  = pair.Midd
		Tail  = pair.Tail
	)
	fmtUtil.Fprintf(
		out,
		"%s\t%s\t%s\t%s\t%d-%d\n",
		left.String("+"),
		Head.String("+"), Midd.String("+"),
		pair.Name, pair.Start, pair.End,
	)
	fmtUtil.Fprintf(
		out,
		"%s\t%s\t%s\t%s\t%d-%d\n",
		right.String("-"),
		Midd.String("-"), Tail.String("-"),
		pair.Name, pair.Start, pair.End,
	)
}
