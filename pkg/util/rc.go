package util

import (
	"sync"
)

// ---------- 查表基础 ----------
var complementTable [256]byte

func init() {
	for i := 0; i < 256; i++ {
		complementTable[i] = byte(i)
	}
	complementTable['A'] = 'T'
	complementTable['T'] = 'A'
	complementTable['C'] = 'G'
	complementTable['G'] = 'C'
	complementTable['a'] = 't'
	complementTable['t'] = 'a'
	complementTable['c'] = 'g'
	complementTable['g'] = 'c'
	complementTable['N'] = 'N'
	complementTable['n'] = 'n'
}

// ---------- 基础版本 ----------
func ReverseComplementBytes(seq []byte) []byte {
	n := len(seq)
	if n == 0 {
		return []byte{}
	}
	out := make([]byte, n)
	for i := range n {
		out[i] = complementTable[seq[n-1-i]]
	}
	return out
}

// ---------- 优化：直接读取 string（避免 []byte 拷贝） ----------
func ReverseComplement2Bytes(s string) []byte {
	n := len(s)
	if n == 0 {
		return []byte{}
	}
	out := make([]byte, n)
	for i := range n {
		out[i] = complementTable[s[n-1-i]]
	}
	return out
}

// ---------- 输入输出均为 string 的优化版 ----------
// 建议基于输入类型使用 ReverseComplement2Bytes 或 ReverseComplementBytes，并在需要输出string时使用 string() 转换
func ReverseComplement(s string) string {
	n := len(s)
	if n == 0 {
		return ""
	}
	b := make([]byte, n)
	for i := range n {
		b[i] = complementTable[s[n-1-i]]
	}
	return string(b)
}

// ---------- 带输出缓冲区的版本 ----------
func ReverseComplementBytesTo(seq []byte, out []byte) []byte {
	n := len(seq)
	if n == 0 {
		return out[:0]
	}
	if out == nil || cap(out) < n {
		out = make([]byte, n)
	} else {
		out = out[:n]
	}
	for i := 0; i < n; i++ {
		out[i] = complementTable[seq[n-1-i]]
	}
	return out
}

// ---------- 使用 sync.Pool 复用切片 ----------
var resultPool = sync.Pool{
	New: func() interface{} { return make([]byte, 0, 1024) },
}

func ReverseComplementBytesPool(seq []byte) []byte {
	n := len(seq)
	if n == 0 {
		return []byte{}
	}
	buf := resultPool.Get().([]byte)
	if cap(buf) < n {
		buf = make([]byte, n)
	}
	buf = buf[:n]
	for i := 0; i < n; i++ {
		buf[i] = complementTable[seq[n-1-i]]
	}
	return buf
}

func ReleasePoolResult(buf []byte) {
	resultPool.Put(buf)
}
