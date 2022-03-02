template <typename T> inline std::pair<int, int> split_vector(std::vector<T> &vec, int comm_size, int rank)
{
    int nblock = vec.size() / comm_size;
    int offset = nblock * rank;
    if (rank == comm_size - 1)
        nblock = vec.size() - offset;

    return std::make_pair(offset, nblock);
}
