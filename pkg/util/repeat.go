package util

import (
	"fmt"
	"sort"
)

// FindRepeats 寻找序列中所有长度 >= minLen 的极大重复区域，
// 返回去除了被更长重复完全覆盖的 Feature 列表。
func FindRepeats(seq string, minLen int, name string) []*Feature {
	n := len(seq)
	if n < minLen {
		return nil
	}

	sa := buildSuffixArray(seq)

	rank := make([]int, n)
	for i := range n {
		rank[sa[i]] = i
	}
	lcp := make([]int, n)
	h := 0
	for i := range n {
		r := rank[i]
		if r == 0 {
			continue
		}
		j := sa[r-1]
		for i+h < n && j+h < n && seq[i+h] == seq[j+h] {
			h++
		}
		lcp[r] = h
		if h > 0 {
			h--
		}
	}

	regions := extractMaximalRepeats(sa, lcp, seq, minLen)

	// 转换为 Feature：为每个出现位置创建一条记录
	features := make([]*Feature, 0)
	for _, reg := range regions {
		firstStart := sa[reg.saStart]
		subSeq := seq[firstStart : firstStart+reg.len]
		count := reg.saEnd - reg.saStart
		for k := reg.saStart; k < reg.saEnd; k++ {
			start := sa[k]
			features = append(features, &Feature{
				Chr:   name,
				Start: start,
				End:   start + reg.len,
				Name:  fmt.Sprintf("Repeat:%s:%d", subSeq, count),
				Seq:   subSeq,
			})
		}
	}

	return removeCovered(features)
}

// ---------- 辅助结构与函数 ----------

type repeatRegion struct {
	saStart int // 在 SA 中的起始索引（包含）
	saEnd   int // 在 SA 中的结束索引（不包含）
	len     int // 公共前缀长度
}

// extractMaximalRepeats 利用单调栈找出所有极大重复区间
func extractMaximalRepeats(sa []int, lcp []int, seq string, minLen int) []repeatRegion {
	n := len(lcp)
	// 哨兵，保证栈被清空
	lcp = append(lcp, 0)
	stack := make([]int, 0, n)
	var regions []repeatRegion

	for i := 0; i <= n; i++ {
		for len(stack) > 0 && lcp[stack[len(stack)-1]] > lcp[i] {
			top := stack[len(stack)-1]
			stack = stack[:len(stack)-1]
			h := lcp[top]
			left := 0
			if len(stack) > 0 {
				left = stack[len(stack)-1] + 1
			}
			right := i

			// 出现次数至少 2 次，且长度达到要求
			if h >= minLen && right-left >= 2 {
				// 检查左右字符以判断是否极大
				allLeftSame, allRightSame := true, true
				var leftChar, rightChar byte = 0, 0
				for k := left; k < right; k++ {
					pos := sa[k]
					// 左边界检查
					if pos == 0 {
						allLeftSame = false
					} else {
						ch := seq[pos-1]
						if k == left {
							leftChar = ch
						} else if ch != leftChar {
							allLeftSame = false
						}
					}
					// 右边界检查
					if pos+h >= len(seq) {
						allRightSame = false
					} else {
						ch := seq[pos+h]
						if k == left {
							rightChar = ch
						} else if ch != rightChar {
							allRightSame = false
						}
					}
				}
				// 只要有一侧字符不全相同，就是极大重复
				if !allLeftSame && !allRightSame {
					regions = append(regions, repeatRegion{
						saStart: left,
						saEnd:   right,
						len:     h,
					})
				} else if !allLeftSame && allRightSame {
					// 右侧全相同 -> 可向右扩展，不是极大（除非右侧已到边界）
				} else if allLeftSame && !allRightSame {
					// 左侧全相同 -> 可向左扩展，不是极大
				}
				// 如果两侧都全相同，显然不是极大（可双向扩展）
			}
		}
		if i < n {
			// 维护栈内 LCP 单调递增
			if len(stack) == 0 || lcp[i] >= lcp[stack[len(stack)-1]] {
				stack = append(stack, i)
			} else {
				// 此情况不会发生，因为 lcp[i] >= lcp[stack[len(stack)-1]] 已保证入栈条件
				stack = append(stack, i)
			}
		}
	}
	return regions
}

// removeCovered 实现与 RepeatDeleteCovered 相同的逻辑，去除被完全覆盖且长度更短的重复
func removeCovered(feats []*Feature) []*Feature {
	if len(feats) == 0 {
		return feats
	}
	// 按 Start 升序，Start 相同时 End 升序
	sort.Slice(feats, func(i, j int) bool {
		if feats[i].Start == feats[j].Start {
			return feats[i].End < feats[j].End
		}
		return feats[i].Start < feats[j].Start
	})

	newList := make([]*Feature, 0, len(feats))
	for i, feat1 := range feats {
		length1 := feat1.End - feat1.Start
		keep := true
		for j, feat2 := range feats {
			if i == j {
				continue
			}
			// 如果 feat2 完全覆盖 feat1 且严格更长，则丢弃 feat1
			if feat2.Start <= feat1.Start && feat2.End >= feat1.End && feat2.End-feat2.Start > length1 {
				keep = false
				break
			}
		}
		if keep {
			newList = append(newList, feat1)
		}
	}
	return newList
}

// buildSuffixArray 倍增法构建后缀数组
func buildSuffixArray(s string) []int {
	n := len(s)
	sa := make([]int, n)
	rank := make([]int, n)
	tmp := make([]int, n)

	for i := 0; i < n; i++ {
		sa[i] = i
		rank[i] = int(s[i])
	}

	k := 1
	// 比较器：先比当前rank，再比后k个字符的rank
	cmp := func(i, j int) bool {
		if rank[i] != rank[j] {
			return rank[i] < rank[j]
		}
		ri, rj := -1, -1
		if i+k < n {
			ri = rank[i+k]
		}
		if j+k < n {
			rj = rank[j+k]
		}
		return ri < rj
	}

	for k < n {
		sort.Slice(sa, cmp)
		tmp[sa[0]] = 0
		for i := 1; i < n; i++ {
			if cmp(sa[i-1], sa[i]) {
				tmp[sa[i]] = tmp[sa[i-1]] + 1
			} else {
				tmp[sa[i]] = tmp[sa[i-1]]
			}
		}
		rank, tmp = tmp, rank
		if rank[sa[n-1]] == n-1 {
			break
		}
		k <<= 1
	}
	return sa
}
