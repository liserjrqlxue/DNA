package util

import (
	"fmt"
	"sort"
)

type Feature struct {
	chr   string
	start int
	end   int
	name  string
	seq   string
}

func (f *Feature) String() string {
	return fmt.Sprintf("%s\t%d\t%d\t%s", f.chr, f.start, f.end, f.name)
}

func GetRepeat(seq string, length int) []*Feature {
	var repeat []*Feature
	for i := 0; i < len(seq)-length; i++ {
		targetSeq := seq[i : i+length]
		for j := i + 1; j < len(seq)-length; j++ {
			if targetSeq == seq[j:j+length] || targetSeq == ReverseComplement(seq[j:j+length]) {
				repeat = append(
					repeat,
					&Feature{
						chr:   Name,
						start: i,
						end:   i + length,
					},
					&Feature{
						chr:   Name,
						start: j,
						end:   j + length,
					},
				)
			}
		}
	}
	return repeat
}

// mergeIntervals 合并有交集的区间
func MergeIntervals(intervals []*Feature) []*Feature {
	if len(intervals) == 0 {
		return intervals
	}

	// 按起点排序
	sort.Slice(intervals, func(i, j int) bool {
		return intervals[i].start < intervals[j].start
	})

	merged := make([]*Feature, 0)
	current := intervals[0]

	for _, interval := range intervals[1:] {
		// fmt.Printf("\t%d-%d\n", interval[0], interval[1])
		if interval.start <= current.end { // 有交集
			// 合并区间
			if interval.end > current.end {
				current.end = interval.end
			}
		} else {
			// 没有交集，保存当前区间并更新
			merged = append(merged, current)
			current = interval
		}
	}

	// 添加最后一个区间
	merged = append(merged, current)

	return merged
}

func SumLength(intervals []*Feature) int {
	var sum int
	for _, interval := range intervals {
		sum += interval.end - interval.start
	}
	return sum
}
