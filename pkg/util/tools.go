package util

import (
	"fmt"
	"log"
	"strings"

	"github.com/liserjrqlxue/goUtil/textUtil"
)

/*
func NearestNeighbor() {
	// Tm=deltaH/(A+deltaS+Rln(C/4))-273.15+16.6log[Na+]
	var (
		A   = -0.0108   // kcal*K-1*mol-1 // 常数，为-0.0108 kcal K-1 ᐧ mol-1（退火/解链过程中的螺旋起始）
		R   = 0.00199   // kcal*K-1*mol-1 // 气体常数，为0.00199 kcal K-1 ᐧ mol-1（针对温度按比例调节能量的常数）
		C   = 0.0000005 // M or mol*L-1 // 寡核苷酸浓度，单位为M或mol L-1（我们使用0.0000005，即0.5 μM）
		K2C = -273.15   // 转换因子，用于将预期的K氏温度转换为°C
		Na  = 0.05      // M // 钠离子浓度，单位为M或mol L-1（我们使用0.05，即50 mM）
	)

}
*/

func LoadInputSeq(path string) string {
	var sequence string
	for _, line := range textUtil.File2Array(path) {
		sequence += strings.TrimSpace(line)
	}
	return strings.ToUpper(sequence)
}

func CheckLocalRepeat(region []*Feature, start, end int) bool {
	for _, feature := range region {
		// 1 left up steam
		if start >= feature.end {
			continue
		} else { // start<feature.end
			if start >= feature.start {
				if end > feature.end { // 2 left overlap
					if feature.end-start > outRepeat { // overlap
						return false
					}
					if end-feature.end < outRepeat { // flank
						return false
					}
				} else { // end<=feature.end // 3 cover
					return false
				}
			} else { // start < feature.start
				if end > feature.end { // 4 inner
					if feature.start-start < outRepeat { // flank
						return false
					}
					if end-feature.end < outRepeat { // flank
						return false
					}
				} else { // end<=feature.end
					if feature.start < end { // 5 right overlap
						if end-feature.start > outRepeat { // overlap
							return false
						}
						if feature.start-start < outRepeat { // flank
							return false
						}
					} else { // end<=feature.start // 6 right
						continue
					}
				}
			}
		}
	}
	return true
}

// CheckGC70 检查序列的GC含量，70%以上返回true，否则返回false
func CheckGC70(seq string) bool {
	var gc = 0
	for _, c := range seq {
		if c == 'C' || c == 'G' {
			gc++
		}
	}
	return float64(gc)/float64(len(seq)) >= 0.7
}

func MergeFeatures(seq string, region []*Feature) []*Feature {
	if region == nil {
		return nil
	}
	var (
		newRegion []*Feature
		reg       = &Feature{
			chr:   region[0].chr,
			start: region[0].start,
			end:   region[0].end,
			name:  region[0].name,
		}
		count = 0
	)
	for i := 0; i < len(region)-1; i++ {
		if reg.end > region[i].start {
			reg.end = Max(reg.end, region[i].end)
			count++
		} else {
			reg.name = fmt.Sprintf("Merged:%s:%d", seq[reg.start:reg.end], count)
			newRegion = append(newRegion, reg)
			reg = &Feature{
				chr:   region[i].chr,
				start: region[i].start,
				end:   region[i].end,
				name:  region[i].name,
			}
			count = 1
		}
	}
	newRegion = append(newRegion, reg)
	return newRegion
}

func NotAllGC(seq string) bool {
	for _, c := range seq {
		switch c {
		case 'A', 'T':
			return true
		}
	}
	return false
}

func Max(x, y int) int {
	if y > x {
		return y
	}
	return x
}

func Min(x, y int) int {
	if y < x {
		return y
	}
	return x
}

// from https://forum.golangbridge.org/t/easy-way-for-letter-substitution-reverse-complementary-dna-sequence/20101
// from https://go.dev/play/p/IXI6PY7XUXN
var dnaComplement = strings.NewReplacer(
	"A", "T",
	"T", "A",
	"G", "C",
	"C", "G",
	"a", "t",
	"t", "a",
	"g", "c",
	"c", "g",
)

func Complement(s string) string {
	return dnaComplement.Replace(s)
}

// Reverse returns its argument string reversed rune-wise left to right.
// from https://github.com/golang/example/blob/master/stringutil/reverse.go
func Reverse(r []byte) []byte {
	for i, j := 0, len(r)-1; i < len(r)/2; i, j = i+1, j-1 {
		r[i], r[j] = r[j], r[i]
	}
	return r
}

// ReverseComplement computes the reverse complement of a DNA sequence.
//
// It takes a string representing a DNA sequence as input.
// It returns a string representing the reverse complement of the input sequence.
func ReverseComplement(s string) string {
	return Complement(string(Reverse([]byte(s))))
}

func FindMin(x []float64) (min float64, index int) {
	min = x[0]
	index = 0
	for i, f := range x {
		if min > f {
			min = f
			index = i
		}
	}
	return
}

func InitLog(seq *Seq) {
	log.Printf("Input Seq length:\t%d", seq.Length)
	log.Printf("Segment length:\t%d [%d-%d]", SegLength, SegLengthRange[0], SegLengthRange[1])
	log.Printf("Overlape length:\t%d [%d-%d]", PrimerLength, PrimerRange[0], PrimerRange[1])
	log.Printf("Tm Range:\t\t%.1f [%.1f,%.1f]", TmRange[2], TmRange[0], TmRange[1])
	log.Printf(
		"Expect primerPair:\t%.1f [%.1f,%.1f]",
		float64(seq.Length-PrimerLength)/float64(SegPairLength-PrimerLength),
		float64(seq.Length-PrimerRange[0])/float64(SegLengthRange[1]-PrimerRange[0])/2,
		float64(seq.Length-PrimerRange[1])/float64(SegLengthRange[0]-PrimerRange[1])/2,
	)
}
