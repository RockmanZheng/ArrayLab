using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

namespace ArrayLabCore
{
    public class LArray<T>
    {
        private T[] _data;
        private T _scalar;
        private int _rank;
        public int dimension { get { return _rank; } }
        private int[] _sizes;
        public int[] shape { get { return _sizes; } }
        // multiplicative cumulation of _sizes
        // used for indexing
        private int[] _mul_sizes;

        // internal helper functions
        private void _set_mul_sizes()
        {
            _mul_sizes = new int[_sizes.Length];
            _mul_sizes[_mul_sizes.Length - 1] = 1;
            for (var i = _mul_sizes.Length - 2; i >= 0; i--)
            {
                _mul_sizes[i] = _sizes[i + 1] * _mul_sizes[i + 1];
            }
        }
        private int _get_index(params int[] idxs)
        {
            if (idxs.Length != _rank)
            {
                throw new ArgumentException("Number of indices mismatches with array rank.");
            }
            for (var i = 0; i < idxs.Length; i++)
            {
                // deal with negative idx
                idxs[i] = _trim_index(idxs[i], i);
                // check out of range idx
                if (idxs[i] < 0 || idxs[i] > _sizes[i])
                {
                    throw new IndexOutOfRangeException();
                }
            }
            int idx = 0;
            for (var i = 0; i < idxs.Length; i++)
            {
                idx += idxs[i] * _mul_sizes[i];
            }
            return idx;
        }
        private int _trim_index(int k, int dim)
        {
            var idx = k < 0 ? k % _sizes[dim] + _sizes[dim] : k;
            if (idx > _sizes[dim])
            {
                throw new IndexOutOfRangeException();
            }
            return idx;
        }
        private int _trim_index(int k)
        {
            var idx = k < 0 ? k % _data.Length + _data.Length : k;
            if (idx > _data.Length)
            {
                throw new IndexOutOfRangeException();
            }
            return idx;
        }
        private Tuple<int, int> _parse_slice_expression(string expr)
        {
            if (!expr.Contains(':'))
            {
                throw new ArgumentException("Not a slice expression.");
            }
            var tokens = expr.Split(':');
            if (tokens.Length > 2)
            {
                throw new ArgumentException("Invalid slice expression. Too many ':'");
            }
            int start;
            int end;
            if (tokens[0] == "")
            {
                start = 0;
            }
            else
            {
                start = _trim_index(int.Parse(tokens[0]));
            }
            if (tokens[1] == "")
            {
                end = _data.Length;
            }
            else
            {
                end = _trim_index(int.Parse(tokens[1]));
            }
            if (start == _data.Length)
            {
                throw new IndexOutOfRangeException();
            }
            if (start == end)
            {
                throw new ArgumentException("Invalid slice expression: collapsed dimension.");
            }
            return Tuple.Create(start, end);
        }
        private Tuple<int, int> _parse_slice_expression(string expr, int dim)
        {
            if (!expr.Contains(':'))
            {
                throw new ArgumentException("Not a slice expression.");
            }
            var tokens = expr.Split(':');
            if (tokens.Length > 2)
            {
                throw new ArgumentException("Invalid slice expression. Too many ':'");
            }
            int start;
            int end;
            if (tokens[0] == "")
            {
                start = 0;
            }
            else
            {
                start = _trim_index(int.Parse(tokens[0]), dim);
            }
            if (tokens[1] == "")
            {
                end = _sizes[dim];
            }
            else
            {
                end = _trim_index(int.Parse(tokens[1]), dim);
            }
            if (start == _sizes[dim])
            {
                throw new IndexOutOfRangeException();
            }
            if (start == end)
            {
                throw new ArgumentException("Invalid slice expression: collapsed dimension.");
            }
            return Tuple.Create(start, end);
        }

        // constructors
        public LArray(T data)
        {
            _rank = 0;
            _data = new T[0];
            _sizes = new int[0];
            _scalar = data;
            _mul_sizes = new int[0];
        }
        public LArray(Array data)
        {
            _rank = data.Rank;
            var len = 1;
            _sizes = new int[_rank];
            for (var i = 0; i < _rank; i++)
            {
                _sizes[i] = data.GetLength(i);
                len *= _sizes[i];
            }
            _set_mul_sizes();
            _data = new T[len];
            //Array.Copy(data.Cast<T>().ToArray(), _data, len);
            Buffer.BlockCopy(data, 0, _data, 0, len * Marshal.SizeOf(typeof(T)));
        }
        public LArray(params int[] dims)
        {
            for (var i = 0; i < dims.Length; i++)
            {
                if (dims[i] <= 0)
                {
                    throw new ArgumentException("Dimension must be positive number.");
                }
            }
            _rank = dims.Length;
            //dims.CopyTo(_sizes,0);
            if (_rank == 0)
            {
                _data = new T[0];
                _sizes = new int[0];
                //_scalar = 0;
                _mul_sizes = new int[0];
            }
            else
            {
                var len = 1;
                _sizes = new int[_rank];
                for (var i = 0; i < _rank; i++)
                {
                    _sizes[i] = dims[i];
                    len *= _sizes[i];
                }
                _set_mul_sizes();
                _data = new T[len];
            }
        }
        public LArray(T val, params int[] dims)
        {
            for (var i = 0; i < dims.Length; i++)
            {
                if (dims[i] <= 0)
                {
                    throw new ArgumentException("Dimension must be positive number.");
                }
            }
            _rank = dims.Length;
            //dims.CopyTo(_sizes,0);
            if (_rank == 0)
            {
                _data = new T[0];
                _sizes = new int[0];
                //_scalar = 0;
                _mul_sizes = new int[0];
            }
            else
            {
                var len = 1;
                _sizes = new int[_rank];
                for (var i = 0; i < _rank; i++)
                {
                    _sizes[i] = dims[i];
                    len *= _sizes[i];
                }
                _set_mul_sizes();
                _data = new T[len];
            }
            // Assignment
            for (var i = 0; i < _data.Length; i++)
            {
                _data[i] = val;
            }
        }
        public LArray(LArray<T> other)
        {
            _rank = other._rank;
            _sizes = new int[other._sizes.Length];
            other._sizes.CopyTo(_sizes, 0);
            _data = new T[other._data.Length];
            other._data.CopyTo(_data, 0);
            _scalar = other._scalar;
        }

        // indexers
        public T this[int idx]
        {
            get
            {
                idx = _trim_index(idx);
                return _data[idx];
            }
            set
            {
                idx = _trim_index(idx);
                _data[idx] = value;
            }
        }
        public LArray<T> this[string index]
        {
            get
            {
                var tokens = index.Split(',');
                // slice support
                if (tokens.Length > _rank)
                {
                    throw new ArgumentException("Index out of dimension.");
                }
                if (tokens.Length != _rank && tokens.Length != 1)
                {
                    throw new ArgumentException("Index dimensions mismatch with data dimensions");
                }
                if (tokens.Length == 1)
                {
                    var token = tokens[0];
                    // slice expression
                    if (token.Contains(':'))
                    {
                        var tmp = _parse_slice_expression(token);
                        var start = tmp.Item1;
                        var end = tmp.Item2;
                        if (end == start)
                        {
                            return new LArray<T>();
                        }
                        var slice = new LArray<T>(Math.Abs(end - start));
                        if (end > start)
                        {
                            Array.Copy(_data, start, slice._data, 0, end - start);
                        }
                        else
                        {
                            Array.Copy(_data, end, slice._data, 0, start - end);
                            slice._data = slice._data.Reverse().ToArray();
                        }
                        return slice;
                    }
                    else// single index
                    {
                        var idx = _trim_index(int.Parse(token));
                        return new LArray<T>(_data[idx]);
                    }
                }
                else// multi-dimensional array situation
                {
                    int rank = tokens.Length;
                    var starts = new List<int>(rank);
                    var ends = new List<int>(rank);
                    var dims = new List<int>(rank);

                    for (var i = 0; i < rank; i++)
                    {
                        var token = tokens[i];
                        // slice expression
                        if (token.Contains(':'))
                        {
                            var tmp = _parse_slice_expression(token, i);
                            starts.Add(tmp.Item1);
                            ends.Add(tmp.Item2);
                            dims.Add(Math.Abs(ends[i] - starts[i]));
                        }
                        else// single index
                        {
                            var idx = _trim_index(int.Parse(token), i);
                            if (idx == _sizes[i])
                            {
                                throw new IndexOutOfRangeException();
                            }
                            starts.Add(idx);
                            ends.Add(idx + 1);
                            dims.Add(1);
                        }
                    }

                    rank = dims.Count;
                    var subarray = new LArray<T>(dims.ToArray());
                    // idxs[k] stores the current working index of k-th dimension of _data
                    var idxs = new int[rank];
                    var subarray_idxs = new int[rank];
                    // initialize idxs
                    starts.CopyTo(0, idxs, 0, rank);
                    int len = subarray._sizes[rank - 1];
                    var buf = new T[len];
                    while (true)
                    {
                        // check if we are done copying slices
                        if (idxs[0] == ends[0]) break;
                        int start = _get_index(idxs);

                        int subarray_start = subarray._get_index(subarray_idxs);
                        // copy slice
                        if (starts[rank - 1] < ends[rank - 1])
                        {
                            Array.Copy(_data, start, subarray._data, subarray_start, len);
                        }
                        else
                        {
                            // reverse order
                            Array.Copy(_data, start - len + 1, buf, 0, len);
                            Array.Copy(buf.Reverse().ToArray(), 0, subarray._data, subarray_start, len);
                        }
                        // increment (or decrement)
                        if (starts[rank - 2] < ends[rank - 2]) idxs[rank - 2]++;
                        else idxs[rank - 2]--;
                        subarray_idxs[rank - 2]++;

                        // make sure all indices are updated
                        for (var i = rank - 2; i > 0; i--)
                        {
                            // if we are done looping i-th dimension
                            if (idxs[i] == ends[i])
                            {
                                idxs[i] = starts[i];
                                if (starts[i - 1] < ends[i - 1]) idxs[i - 1]++;
                                else idxs[i - 1]--;
                            }
                            if (subarray_idxs[i] == subarray._sizes[i])
                            {
                                subarray_idxs[i] = 0;
                                subarray_idxs[i - 1]++;
                            }
                        }
                    }
                    return subarray;
                }
            }
            set
            {
                var tokens = index.Split(',');
                // slice support
                if (tokens.Length > _rank)
                {
                    throw new ArgumentException("Index out of dimension.");
                }
                if (tokens.Length != _rank && tokens.Length != 1)
                {
                    throw new ArgumentException("Index dimensions mismatch with data dimensions");
                }
                // flat assignment
                if (tokens.Length == 1)
                {
                    var token = tokens[0];
                    // slice expression
                    if (token.Contains(':'))
                    {
                        var tmp = _parse_slice_expression(token);
                        var start = tmp.Item1;
                        var end = tmp.Item2;
                        if (value._data.Length != Math.Abs(end - start))
                        {
                            throw new ArgumentException("Assigned array length mismatches with slice length.");
                        }
                        if (start == end)
                        {
                            throw new ArgumentException("Empty sclice, useless assignment.");
                        }
                        if (start < end)
                        {
                            Array.Copy(value._data, 0, _data, start, end - start);
                        }
                        else
                        {
                            Array.Copy(value._data.Reverse().ToArray(), 0, _data, end, start - end);
                        }
                    }
                    else// single index
                    {
                        // verify value is indeed scalar
                        if (value._rank != 0)
                        {
                            throw new ArgumentException("Value is not scalar but index indicates scalar assignment.");
                        }
                        var idx = _trim_index(int.Parse(token));
                        _data[idx] = value._scalar;
                    }
                }
                else// sub-array support
                {
                    int rank = tokens.Length;
                    var starts = new List<int>(rank);
                    var ends = new List<int>(rank);
                    var dims = new List<int>(rank);

                    for (var i = 0; i < rank; i++)
                    {
                        var token = tokens[i];
                        // slice expression
                        if (token.Contains(':'))
                        {
                            var tmp = _parse_slice_expression(token, i);
                            starts.Add(tmp.Item1);
                            ends.Add(tmp.Item2);
                            dims.Add(Math.Abs(ends[i] - starts[i]));
                        }
                        else// single index
                        {
                            var idx = _trim_index(int.Parse(token), i);
                            if (idx == _sizes[i])
                            {
                                throw new IndexOutOfRangeException();
                            }
                            starts.Add(idx);
                            ends.Add(idx + 1);
                            dims.Add(1);
                        }
                    }
                    // Verify that assigned array has the same shape with sub-array
                    for (var i = 0; i < rank; i++)
                    {
                        if (dims[i] != value._sizes[i])
                        {
                            throw new ArgumentException("Assigned array mismatches sub-array in shape.");
                        }
                    }
                    // idxs[k] stores the current working index of k-th dimension of _data
                    var idxs = new int[rank];
                    var subarray_idxs = new int[rank];
                    // initialize idxs
                    starts.CopyTo(0, idxs, 0, rank);
                    int len = value._sizes[rank - 1];
                    var buf = new T[len];
                    while (true)
                    {
                        // check if we are done copying slices
                        if (idxs[0] == ends[0]) break;

                        int start = _get_index(idxs);
                        int subarray_start = value._get_index(subarray_idxs);
                        // copy slice
                        if (starts[rank - 1] < ends[rank - 1])
                        {
                            Array.Copy(value._data, subarray_start, _data, start, len);
                        }
                        else
                        {
                            // reverse order
                            Array.Copy(value._data, subarray_start, buf, 0, len);
                            Array.Copy(buf.Reverse().ToArray(), 0, _data, start - len + 1, len);
                        }
                        // increment (or decrement)
                        if (starts[rank - 2] < ends[rank - 2]) idxs[rank - 2]++;
                        else idxs[rank - 2]--;
                        subarray_idxs[rank - 2]++;

                        // make sure all indices are updated
                        for (var i = rank - 2; i > 0; i--)
                        {
                            // if we are done looping i-th dimension
                            if (idxs[i] == ends[i])
                            {
                                idxs[i] = starts[i];
                                if (starts[i - 1] < ends[i - 1]) idxs[i - 1]++;
                                else idxs[i - 1]--;
                            }
                            if (subarray_idxs[i] == value._sizes[i])
                            {
                                subarray_idxs[i] = 0;
                                subarray_idxs[i - 1]++;
                            }
                        }
                    }
                }
            }
        }

        public override string ToString()
        {
            string result = "";
            // scalar situation
            if (_rank == 0)
            {
                result += _scalar.ToString();
                return result;
            }

            var idxs = new int[_rank];
            var count = _rank;
            while (true)
            {
                for (var i = 0; i < count; i++)
                {
                    result += "[";
                }
                count = 0;
                var idx = _get_index(idxs);
                result += _data[idx].ToString();
                // increment
                idxs[_rank - 1]++;
                for (var i = _rank - 1; i > 0; i--)
                {
                    if (idxs[i] == _sizes[i])
                    {
                        idxs[i] = 0;
                        idxs[i - 1]++;
                        result += "]";
                        count++;
                    }
                }

                if (idxs[0] == _sizes[0])
                {
                    result += "]";
                    break;
                }
                result += ",";
            }
            return result;
        }

        // equality operators
        public static bool operator ==(LArray<T> a, LArray<T> b)
        {
            if (a._rank != b._rank) return false;
            for (var i = 0; i < a._rank; i++)
            {
                if (a._sizes[i] != b._sizes[i]) return false;
            }
            if (a._rank == 0)
            {
                return EqualityComparer<T>.Default.Equals(a._scalar, b._scalar);
            }
            for (var i = 0; i < a._data.Length; i++)
            {
                if (!EqualityComparer<T>.Default.Equals(a[i], b[i])) return false;
            }
            return true;
        }
        public static bool operator !=(LArray<T> a, LArray<T> b)
        {
            if (a._rank != b._rank) return true;
            for (var i = 0; i < a._rank; i++)
            {
                if (a._sizes[i] != b._sizes[i]) return true;
            }
            if (a._rank == 0)
            {
                return !EqualityComparer<T>.Default.Equals(a._scalar, b._scalar);
            }
            for (var i = 0; i < a._data.Length; i++)
            {
                if (!EqualityComparer<T>.Default.Equals(a[i], b[i])) return true;
            }
            return false;
        }
        public static bool operator ==(LArray<T> a, T b)
        {
            if (a._rank != 0) return false;
            return EqualityComparer<T>.Default.Equals(a._scalar, b);
        }
        public static bool operator !=(LArray<T> a, T b)
        {
            return a._rank != 0 || !EqualityComparer<T>.Default.Equals(a._scalar, b);
        }
        public static bool operator ==(T a, LArray<T> b)
        {
            return b._rank == 0 && EqualityComparer<T>.Default.Equals(a, b._scalar);
        }
        public static bool operator !=(T a, LArray<T> b)
        {
            return b._rank != 0 || !EqualityComparer<T>.Default.Equals(a, b._scalar);
        }

        // to linq multi-dimensional array
        public Array ToArray()
        {
            if (_rank == 0)
            {
                throw new ArgumentException("Cannot convert scalar to array");
            }
            var a = Array.CreateInstance(typeof(T), _sizes);
            Buffer.BlockCopy(_data, 0, a, 0, Marshal.SizeOf(typeof(T)) * _data.Length);
            return a;
        }

        public LArray<T> reshape(params int[] dims)
        {
            var len = 1;
            for (var i = 0; i < dims.Length; i++)
            {
                len *= dims[i];
            }
            if (len != _data.Length)
            {
                throw new ArgumentException("Cannot reshape array due to mismatching length");
            }
            var new_array = new LArray<T>(dims);
            _data.CopyTo(new_array._data, 0);
            return new_array;
        }

        // reverse dimension order
        public LArray<T> reverse(int dim)
        {
            if (dim < 0 || dim >= _rank)
            {
                throw new ArgumentException("Cannot reverse array due to invalid dimension.");
            }
            var result = new LArray<T>(_sizes);
            var idxs = new int[_rank];
            //var result_idxs = new int[_rank];
            if (dim < _rank - 1) idxs[dim] = _sizes[dim] - 1;
            int end_0;
            if (dim == 0) end_0 = -1;
            else end_0 = _sizes[0];
            var result_idx = 0;
            var len = _sizes[_rank - 1];
            while (true)
            {
                if (idxs[0] == end_0) break;
                int idx;
                var buf = new T[len];
                idx = _get_index(idxs);
                if (dim == _rank - 1)
                {
                    Array.Copy(_data, idx, buf, 0, len);
                    Array.Copy(buf.Reverse().ToArray(), 0, result._data, result_idx, len);
                }
                else
                {
                    Array.Copy(_data, idx, result._data, result_idx, len);
                }
                // update index
                result_idx += len;
                if (dim == _rank - 2) idxs[_rank - 2]--;
                else idxs[_rank - 2]++;
                for (var i = _rank - 2; i > 0; i--)
                {
                    if (dim == i)
                    {
                        if (idxs[i] < 0)
                        {
                            idxs[i] = _sizes[i] - 1;
                            idxs[i - 1]++;
                        }
                    }
                    else
                    {
                        if (idxs[i] >= _sizes[i])
                        {
                            idxs[i] = 0;
                            if (i - 1 == dim) idxs[i - 1]--;
                            else idxs[i - 1]++;
                        }
                    }
                }
            }
            return result;
        }

        // basic operators
        public static LArray<T> operator +(LArray<T> a, LArray<T> b)
        {
            if (!a._sizes.SequenceEqual(b._sizes))
            {
                throw new ArgumentException("Cannot add 2 arrays with different shapes");
            }
            var c = new LArray<T>(a._sizes);
            for (var i = 0; i < c._data.Length; i++)
            {
                c[i] = (dynamic)a[i] + b[i];
            }
            return c;
        }
        public static LArray<T> operator +(LArray<T> a, T b)
        {
            var c = new LArray<T>(a._sizes);
            for (var i = 0; i < c._data.Length; i++)
            {
                c[i] = (dynamic)a[i] + b;
            }
            return c;
        }
        public static LArray<T> operator +(T b, LArray<T> a)
        {
            return a + b;
        }
        public static LArray<T> operator -(LArray<T> a, LArray<T> b)
        {
            if (!a._sizes.SequenceEqual(b._sizes))
            {
                throw new ArgumentException("Cannot subract 2 arrays with different shapes");
            }
            var c = new LArray<T>(a._sizes);
            for (var i = 0; i < c._data.Length; i++)
            {
                c[i] = (dynamic)a[i] - b[i];
            }
            return c;
        }
        public static LArray<T> operator -(LArray<T> a, T b)
        {
            var c = new LArray<T>(a._sizes);
            for (var i = 0; i < c._data.Length; i++)
            {
                c[i] = (dynamic)a[i] - b;
            }
            return c;
        }
        public static LArray<T> operator -(LArray<T> a)
        {
            var c = new LArray<T>(a._sizes);
            for (var i = 0; i < c._data.Length; i++)
            {
                c[i] = -(dynamic)a[i];
            }
            return c;
        }
        public static LArray<T> operator -(T b, LArray<T> a)
        {
            var c = new LArray<T>(a._sizes);
            for (var i = 0; i < c._data.Length; i++)
            {
                c[i] = (dynamic)b - a[i];
            }
            return c;
        }
        public static LArray<T> operator *(LArray<T> a, LArray<T> b)
        {
            if (!a._sizes.SequenceEqual(b._sizes))
            {
                throw new ArgumentException("Cannot multiply 2 arrays with different shapes");
            }
            var c = new LArray<T>(a._sizes);
            for (var i = 0; i < c._data.Length; i++)
            {
                c[i] = (dynamic)a[i] * b[i];
            }
            return c;
        }
        public static LArray<T> operator *(LArray<T> a, T b)
        {
            var c = new LArray<T>(a._sizes);
            for (var i = 0; i < c._data.Length; i++)
            {
                c[i] = (dynamic)a[i] * b;
            }
            return c;
        }
        public static LArray<T> operator *(T b, LArray<T> a)
        {
            return a * b;
        }
        public static LArray<T> operator /(LArray<T> a, LArray<T> b)
        {
            if (!a._sizes.SequenceEqual(b._sizes))
            {
                throw new ArgumentException("Cannot divide 2 arrays with different shapes");
            }
            var c = new LArray<T>(a._sizes);
            for (var i = 0; i < c._data.Length; i++)
            {
                c[i] = (dynamic)a[i] / b[i];
            }
            return c;
        }
        public static LArray<T> operator /(LArray<T> a, T b)
        {
            var c = new LArray<T>(a._sizes);
            for (var i = 0; i < c._data.Length; i++)
            {
                c[i] = (dynamic)a[i] / b;
            }
            return c;
        }
        public static LArray<T> operator /(T b, LArray<T> a)
        {
            var c = new LArray<T>(a._sizes);
            for (var i = 0; i < c._data.Length; i++)
            {
                c[i] = (dynamic)b / a[i];
            }
            return c;
        }

        // conversion operators
        public static explicit operator T(LArray<T> a)
        {
            if (a._rank != 0)
            {
                throw new InvalidCastException("Cannot convert array into scalar");
            }
            return a._scalar;
        }
        public static explicit operator Array(LArray<T> a)
        {
            return a.ToArray();
        }

        // operations
        public LArray<T> sum(int dim)
        {
            if (dim < 0 || dim >= _rank)
            {
                throw new ArgumentException("Cannot sum over due to invalid dimension");
            }
            var dims = new int[_rank];
            _sizes.CopyTo(dims, 0);
            dims[dim] = 1;
            var result = new LArray<T>(dims);
            var index_start = "";
            var index_end = "";
            for (var i = 0; i < dim; i++) index_start += ":,";
            for (var i = dim + 1; i < _rank; i++) index_end += ",:";
            string index;
            for (var i = 0; i < _sizes[dim]; i++)
            {

                index = index_start + i.ToString() + index_end;
                result += this[index];
            }
            return result;
        }
        public LArray<T> cumsum(int dim)
        {
            if (dim < 0 || dim >= _rank)
            {
                throw new ArgumentException("Cannot sum over due to invalid dimension");
            }
            var result = new LArray<T>(_sizes);
            var index_start = "";
            var index_end = "";
            for (var i = 0; i < dim; i++) index_start += ":,";
            for (var i = dim + 1; i < _rank; i++) index_end += ",:";
            string index_1;
            string index_2;
            index_1 = index_start + "0" + index_end;
            result[index_1] = this[index_1];
            for (var i = 1; i < _sizes[dim]; i++)
            {
                index_1 = index_start + i.ToString() + index_end;
                index_2 = index_start + (i - 1).ToString() + index_end;
                result[index_1] = this[index_1] + result[index_2];
            }
            return result;
        }
        public LArray<T> mul(int dim)
        {
            if (dim < 0 || dim >= _rank)
            {
                throw new ArgumentException("Cannot sum over due to invalid dimension");
            }
            var dims = new int[_rank];
            _sizes.CopyTo(dims, 0);
            dims[dim] = 1;
            var index_start = "";
            var index_end = "";
            for (var i = 0; i < dim; i++) index_start += ":,";
            for (var i = dim + 1; i < _rank; i++) index_end += ",:";
            string index;
            index = index_start + "0" + index_end;
            var result = new LArray<T>(this[index]);
            for (var i = 1; i < _sizes[dim]; i++)
            {
                index = index_start + i.ToString() + index_end;
                result *= this[index];
            }
            return result;
        }
        public LArray<T> cummul(int dim)
        {
            if (dim < 0 || dim >= _rank)
            {
                throw new ArgumentException("Cannot sum over due to invalid dimension");
            }
            var result = new LArray<T>(_sizes);
            var index_start = "";
            var index_end = "";
            for (var i = 0; i < dim; i++) index_start += ":,";
            for (var i = dim + 1; i < _rank; i++) index_end += ",:";
            string index_1;
            string index_2;
            index_1 = index_start + "0" + index_end;
            result[index_1] = this[index_1];
            for (var i = 1; i < _sizes[dim]; i++)
            {
                index_1 = index_start + i.ToString() + index_end;
                index_2 = index_start + (i - 1).ToString() + index_end;
                result[index_1] = this[index_1] * result[index_2];
            }
            return result;
        }
    }
}
