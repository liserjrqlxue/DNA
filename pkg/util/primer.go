package util

import (
	"fmt"
	"math"
)

const Grate = 0.8

type Primer struct {
	Name string `json:"name"`
	Seq  string `json:"seq"`

	Length int `json:"length"`
	Start  int `json:"start"`
	End    int `json:"end"`

	ATCG     map[byte]int
	DupLink  map[byte]int
	GC       float64 `json:"gc"`
	Tm       float64 `json:"tm"`
	DupGrate float64
}

// NewPrimer creates a new Primer with the given parameters.
// It takes in the name, sequence, start, and end values to initialize the Primer.
// It returns a pointer to the created Primer.
func NewPrimer(name string, seq string, start, end int) *Primer {
	// Create a new Primer object and assign the values to its fields.
	primer := &Primer{
		Name:   name,
		Seq:    seq[start:end],
		Start:  start,
		End:    end,
		Length: end - start,
	}

	// Calculate additional properties of the Primer.
	primer.Calculate()

	// Return the created Primer.
	return primer
}

// Calculate calculates the ATCG count, duplicate link count,
// GC content, and melting temperature (Tm).
//
// It takes no parameters and returns no value.
func (p *Primer) Calculate() {
	// Initialize the ATCG and DupLink maps
	p.ATCG = make(map[byte]int)
	p.DupLink = make(map[byte]int)

	// Count the occurrences of each nucleotide in the sequence
	for i := 0; i < p.Length; i++ {
		p.ATCG[p.Seq[i]]++
	}

	// Count the occurrences of duplicate links in the sequence
	for i := 1; i < p.Length; i++ {
		if p.Seq[i] == p.Seq[i-1] {
			p.DupLink[p.Seq[i]]++
		}
	}

	// Calculate the GC content
	GC := p.ATCG['G'] + p.ATCG['C']
	p.GC = float64(GC*100) / float64(p.Length)

	// Calculate the melting temperature (Tm)
	p.Tm = p.CalculateTm()
}

const (
	A = 69.3 // 64.9
	B = 41.0
	C = 650.0 // 41*16.4=672.4
	// D = 16.6
	// Na = 0.05
)

// CalculateTm calculates the melting temperature (Tm) based on the GC content and sequence length.
// It takes the GC content and sequence length as parameters and returns the melting temperature (Tm).
// Calculate Tm using the formula: Tm = A + (B * GC/100) - (C / Length)
func (p *Primer) CalculateTm() float64 {
	return A + B*p.GC/100 - C/float64(p.Length)
}

// String returns a string representation of the Primer.
//
// The String function takes a strand parameter and returns a string representation
// of the Primer. If the strand is "-", the sequence is reversed and complemented,
// and the duplication rate is calculated based on the 'C' duplication link. If
// the strand is not "-", the duplication rate is calculated based on the 'G'
// duplication link.
//
// Parameters:
// - strand: a string indicating the strand value
//
// Return type:
// - string: a formatted string representation of the Primer
func (p *Primer) String(strand string) string {
	// log.Printf("%+v\n", p)
	if p == nil {
		return "\t\t\t\t\t\t\t\t\t"
	}
	seq := p.Seq
	if strand == "-" {
		seq = ReverseComplement(seq)
		p.DupGrate = math.Pow(Grate, float64(p.DupLink['C']))
	} else {
		p.DupGrate = math.Pow(Grate, float64(p.DupLink['G']))
	}
	return fmt.Sprintf(
		"%s\t%s\t%d\t%.1f%%\t%.1f℃\t%s\t%d\t%d\t%f\t%d-%d",
		p.Name, seq, p.Length,
		p.GC, p.Tm, strand,
		p.DupLink['G'], p.DupLink['C'], p.DupGrate,
		p.Start, p.End,
	)
}

func (p *Primer) Check(offset int, repeat []*Feature) bool {
	return p.isInRange() && p.isGCZero() && p.CheckLocalRepeat(offset, repeat)
}

func (p *Primer) isGCZero() bool {
	return p.isTailGCZero() && p.isHeadGCZero()
}

func (p *Primer) isInRange() bool {
	//start := TmRange[2]
	start := TmRange[0]
	end := TmRange[1]
	return p.Tm >= start && p.Tm <= end
}

func (p *Primer) isInRangeCandidate() bool {
	start := TmRange[0]
	end := TmRange[1]
	return p.Tm >= start && p.Tm <= end
}

func (p *Primer) isInSubSeqRangeCandidate() bool {
	start := SubSeqTmRange[0]
	end := SubSeqTmRange[1]
	return p.Tm >= start && p.Tm <= end
}

func (p *Primer) isTailGCZero() bool {
	return NotAllGC(p.Seq[p.Length-terminalGC-1:])
}

func (p *Primer) isHeadGCZero() bool {
	return NotAllGC(p.Seq[:terminalGC+1])
}

// CheckLocalRepeat checks if a local repeat is present in the sequence. True for no repeat
func (p *Primer) CheckLocalRepeat(offset int, repeat []*Feature) bool {
	return CheckLocalRepeat(repeat, p.Start+offset, p.End+offset)
}

// 是否unambiguous ACGT
func (p *Primer) IsUnambiguous() bool {
	return ACGT.MatchString(p.Seq)
}
