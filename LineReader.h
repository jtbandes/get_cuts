#if !defined(__cplusplus) || __cplusplus < 201703L
#error "This file requires C++17"
#endif

#include <fstream>
#include <iostream>
#include <vector>

#include "Progress.h"

// Helper class to read a file line by line, and parse values out of the most recently read line.
class LineReader {
    static const size_t MAX_LINE_LENGTH = 1024;

    char* _p = nullptr;  // current position in line
    char* _end = nullptr;  // end of line

    char _buf[MAX_LINE_LENGTH];
    std::ifstream _file;
    Progress _progress; // must be after _file since we use _file during initialization

    void checkEnd() {
        if (_p == _end) {
            throw std::out_of_range("Read past end of line");
        }
    }

public:
    LineReader(const char* filename)
        : _file(filename, std::ios::in | std::ios::ate /* open at end to get file size */)
        , _progress(filename, _file.tellg())
    {
        _file.seekg(0);
        if (!_file) {
            throw std::system_error(errno, std::system_category(), std::string("Error opening ") + filename);
        }
    }

    // Load a new line from the file. Returns true if the operation succeeded, false if the end of the file was reached.
    bool nextLine() {
        _file.getline(_buf, MAX_LINE_LENGTH);
        if (_file.eof()) {
            _progress.finish();
            _p = nullptr;
            _end = nullptr;
            return false;
        }
        else if (_file.fail()) {
            throw std::length_error("Max line length exceeded");
        }
        size_t len = std::strlen(_buf);

        // Update progress with strlen instead of _file.tellg() which is slow because it actually seeks the file :(
        _progress.addBytesRead(len + 1);

        _p = _buf;
        _end = _p + len;
        return bool(_file);
    }

    // True if all characters on the current line have been consumed
    bool usedWholeLine() const {
        return _p == _end;
    }

    // True if the last call to `nextLine()` reached the end of the input file
    bool atEOF() const {
        return _file.eof();
    }

    // Validate that `str` appears next in the line, and consume it
    void skip(const char* str) {
        checkEnd();
        if (std::strcmp(_p, str) != 0) {
            throw std::runtime_error(std::string("Expected ") + str);
        }
        _p += std::strlen(str);
    }

    // Validate that `c` appears next in the line, and consume it
    void skip(char c) {
        checkEnd();
        if (*_p != c) {
            throw std::runtime_error(std::string("Expected ") + c);
        }
        ++_p;
    }

    // Read the next character from the current line without consuming it
    char peek() {
        checkEnd();
        return *_p;
    }

    // Consume and discard the next N whitespace-separated floating-point values
    template<size_t N> void skipDouble() {
        if constexpr (N > 0) {
            readDouble();
            skipDouble<N-1>();
        }
    }

    // Skip whitespace and consume the next floating-point value
    double readDouble() {
        char* end = _end;
        double val = std::strtod(_p, &end);
        if (end == _p) {
            throw std::runtime_error("Unable to read double");
        }
        _p = end;
        return val;
    }

    // Consume and save comma+whitespace-separated floating-point values until the end of the current line
    void readCommaSeparatedDoubles(std::vector<double>* out) {
        while (true) {
            out->push_back(readDouble());
            if (usedWholeLine()) {
                break;
            } else {
                skip(',');
            }
        }
    }
};
