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

	// 提前计算 Feature 总数
	totalFeatures := 0
	for _, reg := range regions {
		totalFeatures += reg.saEnd - reg.saStart
	}

	// 转换为 Feature：为每个出现位置创建一条记录
	features := make([]*Feature, 0, totalFeatures) // 预分配精确容量
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

func removeCovered(feats []*Feature) []*Feature {
	if len(feats) == 0 {
		return feats
	}
	// 排序：Start 升序，Start 相同时 End 降序（长的在前）
	sort.Slice(feats, func(i, j int) bool {
		if feats[i].Start != feats[j].Start {
			return feats[i].Start < feats[j].Start
		}
		return feats[i].End > feats[j].End
	})

	stack := make([]*Feature, 0, len(feats))
	for _, f := range feats {
		for len(stack) > 0 {
			top := stack[len(stack)-1]
			// 栈顶完全覆盖 f（包含相同区间）
			if top.Start <= f.Start && top.End >= f.End {
				f = nil // 标记丢弃
				break
			}
			// f 完全覆盖栈顶
			if f.Start <= top.Start && f.End >= top.End {
				stack = stack[:len(stack)-1] // 弹出栈顶
				continue
			}
			// 互不包含，停止检查栈顶
			break
		}
		if f != nil {
			stack = append(stack, f)
		}
	}
	return stack
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

func FindInvertedRepeats(seq string, minLen int, name string) []*Feature {
	L := len(seq)
	if L < minLen {
		return nil
	}
	rc := reverseComplement(seq)
	// 拼接：S + "#" + RC + "$"
	combined := seq + "#" + rc + "$"
	sa := buildSuffixArray(combined)
	n := len(combined)
	rank := make([]int, n)
	for i := 0; i < n; i++ {
		rank[sa[i]] = i
	}
	lcp := make([]int, n)
	h := 0
	for i := 0; i < n; i++ {
		r := rank[i]
		if r == 0 {
			continue
		}
		j := sa[r-1]
		for i+h < n && j+h < n && combined[i+h] == combined[j+h] {
			h++
		}
		lcp[r] = h
		if h > 0 {
			h--
		}
	}

	// 用单调栈收集跨区重复
	type region struct {
		saStart, saEnd int
		len            int
	}
	var regions []region
	stack := make([]int, 0, n)
	lcpExt := append(lcp, 0) // 哨兵
	for i := 0; i <= n; i++ {
		for len(stack) > 0 && lcpExt[stack[len(stack)-1]] > lcpExt[i] {
			top := stack[len(stack)-1]
			stack = stack[:len(stack)-1]
			h := lcpExt[top]
			left := 0
			if len(stack) > 0 {
				left = stack[len(stack)-1] + 1
			}
			right := i
			if h >= minLen && right-left >= 2 {
				// 检查区间内是否既有 S 区后缀也有 RC 区后缀
				hasS, hasRC := false, false
				for k := left; k < right; k++ {
					pos := sa[k]
					if pos < L {
						hasS = true
					} else if pos > L && pos < 2*L+1 {
						hasRC = true
					}
				}
				if hasS && hasRC {
					// 进一步要求不能全部在同一侧（因为我们要跨区）
					allSameSide := true
					side := -1
					for k := left; k < right; k++ {
						pos := sa[k]
						if pos < L {
							if side == -1 {
								side = 0
							} else if side != 0 {
								allSameSide = false
								break
							}
						} else if pos > L && pos < 2*L+1 {
							if side == -1 {
								side = 1
							} else if side != 1 {
								allSameSide = false
								break
							}
						}
					}
					if !allSameSide { // 真正跨区
						regions = append(regions, region{saStart: left, saEnd: right, len: h})
					}
				}
			}
		}
		if i < n {
			stack = append(stack, i)
		}
	}

	// 生成 Feature：为每个跨区重复的每个出现位置生成 Feature
	var feats []*Feature
	for _, reg := range regions {
		k := reg.len
		for idx := reg.saStart; idx < reg.saEnd; idx++ {
			pos := sa[idx]
			if pos < L { // S 区
				feats = append(feats, &Feature{
					Chr:   name,
					Start: pos,
					End:   pos + k,
					Name:  fmt.Sprintf("InvRepeat:%s:%d", seq[pos:pos+k], reg.saEnd-reg.saStart),
					Seq:   seq[pos : pos+k],
				})
			} else if pos > L && pos < 2*L+1 { // RC 区，映射回正向
				offsetRC := pos - (L + 1)
				newStart := L - offsetRC - k
				feats = append(feats, &Feature{
					Chr:   name,
					Start: newStart,
					End:   newStart + k,
					Name:  fmt.Sprintf("InvRepeat:%s:%d", seq[newStart:newStart+k], reg.saEnd-reg.saStart),
					Seq:   seq[newStart : newStart+k],
				})
			}
		}
	}
	return feats
}

func (s *Seq) CalculatorRepeat8ntWithAll(name string, queries []string) {
	forward := FindRepeats(s.Seq, 8, name)
	inverted := FindInvertedRepeats(s.Seq, 8, name)
	fuzzy := FuzzyMatchKmersWithRC(s.Seq, queries, name)
	all := append(append(forward, inverted...), fuzzy...)
	s.Repeat = removeCovered(all)
}

// FuzzyMatchKmersWithRC 模糊匹配（正向+反向互补）
func FuzzyMatchKmersWithRC(target string, queries []string, name string) []*Feature {
	var features []*Feature
	features = make([]*Feature, 0, len(target)/2)
	const kmerLen = 8

	for _, query := range queries {
		// 正向和反向互补均作为查询
		qList := []string{query, reverseComplement(query)}
		for _, q := range qList {
			if len(q) < kmerLen {
				continue
			}
			for start := 0; start <= len(q)-kmerLen; start++ {
				kmer := q[start : start+kmerLen]
				// 在 target 正链上扫描
				for i := 0; i <= len(target)-kmerLen; i++ {
					// 精确匹配
					if target[i:i+kmerLen] == kmer {
						features = append(features, &Feature{
							Chr: name, Start: i, End: i + kmerLen,
							Name: fmt.Sprintf("Fuzzy:%s:0", kmer),
							Seq:  kmer,
						})
						continue
					}
					// 单错配
					if matchEditDistance1(target[i:i+kmerLen], kmer) {
						features = append(features, &Feature{
							Chr: name, Start: i, End: i + kmerLen,
							Name: fmt.Sprintf("Fuzzy:%s:1", kmer),
							Seq:  target[i : i+kmerLen],
						})
					}
					// 插入缺失
					if i+kmerLen+1 <= len(target) && matchIndel(target[i:i+kmerLen+1], kmer) {
						features = append(features, &Feature{
							Chr: name, Start: i, End: i + kmerLen + 1,
							Name: fmt.Sprintf("Fuzzy:%s:1", kmer),
							Seq:  target[i : i+kmerLen+1],
						})
					}
					if kmerLen > 1 && i+kmerLen-1 <= len(target) && matchIndel(kmer, target[i:i+kmerLen-1]) {
						features = append(features, &Feature{
							Chr: name, Start: i, End: i + kmerLen - 1,
							Name: fmt.Sprintf("Fuzzy:%s:1", kmer),
							Seq:  target[i : i+kmerLen-1],
						})
					}
				}
			}
		}
	}
	return features
}

// reverseComplement 返回 DNA 序列的反向互补
func reverseComplement(s string) string {
	// 需处理大小写、简并碱基等，这里给出简单版
	complement := map[byte]byte{'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
	res := make([]byte, len(s))
	for i := 0; i < len(s); i++ {
		c := s[len(s)-1-i]
		if v, ok := complement[c]; ok {
			res[i] = v
		} else {
			res[i] = c // 非标准碱基保留
		}
	}
	return string(res)
}

// matchEditDistance1 检查 a 和 b 是否恰好相差 1 个错配（长度必须相等）
func matchEditDistance1(a, b string) bool {
	if len(a) != len(b) {
		return false
	}
	mismatches := 0
	for i := 0; i < len(a); i++ {
		if a[i] != b[i] {
			mismatches++
			if mismatches > 1 {
				return false
			}
		}
	}
	return mismatches == 1
}

// matchIndel 检查 longer 是否比 shorter 多一个字符，且其余部分完全匹配。
// 即 longer 删除某一个字符后等于 shorter。
func matchIndel(longer, shorter string) bool {
	if len(longer) != len(shorter)+1 {
		return false
	}
	for i := 0; i < len(longer); i++ {
		// 尝试跳过 longer 的第 i 个字符
		s := longer[:i] + longer[i+1:]
		if s == shorter {
			return true
		}
	}
	return false
}
