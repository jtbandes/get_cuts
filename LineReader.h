#pragma once

#if !defined(__cplusplus) || __cplusplus < 201703L
#error "This file requires C++17"
#endif

#include <cstdio>
#include <cstring>
#include <memory>
#include <system_error>
#include <vector>

#include "Progress.h"

size_t getFileSize(std::FILE* file) {
    if (file) {
        if (std::fseek(file, 0, SEEK_END) != 0) {
            throw std::system_error(errno, std::system_category(), "Error seeking to end");
        }
        return std::ftell(file);
    }
    return 0;
}

// Helper class to read a file line by line, and parse values out of the most recently read line.
class LineReader {
    static const size_t MAX_LINE_LENGTH = 1024;

    char* _p = nullptr;  // current position in line
    char* _end = nullptr;  // end of line

    char _buf[MAX_LINE_LENGTH];
    std::unique_ptr<std::FILE, decltype(&std::fclose)> _file;
    Progress _progress; // must be after _file since we use _file during initialization

    void checkEnd() {
        if (_p == _end) {
            throw std::out_of_range("Read past end of line");
        }
    }

public:
    LineReader(const char* filename)
        : _file(std::fopen(filename, "r"), std::fclose)
        , _progress(filename, getFileSize(_file.get()))
    {
        if (!_file) {
            throw std::system_error(errno, std::system_category(), std::string("Error opening ") + filename);
        }
        std::rewind(_file.get());
    }

    // Load a new line from the file. Returns true if the operation succeeded, false if the end of the file was reached.
    bool nextLine() {
        if (!std::fgets(_buf, MAX_LINE_LENGTH, _file.get())) {
            if (atEOF()) {
                _progress.finish();
                _p = nullptr;
                _end = nullptr;
                return false;
            }
            throw std::system_error(errno, std::system_category(), "Error reading line from file");
        }

        size_t len = std::strlen(_buf);
        _progress.addBytesRead(len);

        if (len == 0) {
            _p = nullptr;
            _end = nullptr;
            return false;
        }
        if (_buf[len - 1] == '\n') {
            _buf[len - 1] = 0;  // allows for later use of strcmp()
            --len;
        } else if (!atEOF()) {
            throw std::length_error("Max line length exceeded");
        }

        _p = _buf;
        _end = _p + len;
        return true;
    }

    // True if all characters on the current line have been consumed
    bool usedWholeLine() const {
        return _p == _end;
    }

    // True if the last call to `nextLine()` reached the end of the input file
    bool atEOF() const {
        return std::feof(_file.get());
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
