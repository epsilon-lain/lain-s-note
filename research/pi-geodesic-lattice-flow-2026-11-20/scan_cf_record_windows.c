/* Keep two-sided windows around large digits in a newline-separated CF stream.
   Usage: scanner THRESHOLD ZERO_BASED_FIRST_INDEX < iteration.txt > windows.jsonl
   Emit the count and maximum on stderr. The beginning and end of the file
   deliberately have shorter windows. Never treat a missing neighbor as zero. */
#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

enum { radius = 64, capacity = 2 * radius + 1 };

static int parse_uint(const char *text, uint64_t *answer) {
    if (*text < '0' || *text > '9') {
        return 0;
    }
    char *end;
    errno = 0;
    uintmax_t value = strtoumax(text, &end, 10);
    if (errno || value > UINT64_MAX) {
        return 0;
    }
    if (*end == '\r') {
        ++end;
    }
    if (*end == '\n') {
        ++end;
    }
    *answer = (uint64_t)value;
    return *end == '\0';
}

static void emit_window(const uint64_t ring[capacity], uint64_t first,
                        uint64_t last, uint64_t index, uint64_t threshold) {
    uint64_t digit = ring[index % capacity];
    if (index == 0 || digit < threshold) {
        return;
    }
    if (first == 0) {
        first = 1;  /* Exclude a_0 from the reverse tail. */
    }
    printf("{\"index\":%" PRIu64 ",\"digit\":%" PRIu64
           ",\"left_nearest_first\":[", index, digit);
    for (uint64_t j = 1; j <= radius && j <= index - first; ++j) {
        printf("%s%" PRIu64, j == 1 ? "" : ",", ring[(index - j) % capacity]);
    }
    printf("],\"right_nearest_first\":[");
    for (uint64_t j = 1; j <= radius && j <= last - index; ++j) {
        printf("%s%" PRIu64, j == 1 ? "" : ",", ring[(index + j) % capacity]);
    }
    puts("]}");
}

int main(int argc, char **argv) {
    uint64_t threshold, first;
    if (argc != 3 || !parse_uint(argv[1], &threshold) || threshold == 0
        || !parse_uint(argv[2], &first)) {
        fputs("Usage: scanner THRESHOLD ZERO_BASED_FIRST_INDEX\n", stderr);
        return 2;
    }
    uint64_t ring[capacity], index = first, maximum = 0, digit;
    char *line = NULL;
    size_t allocated = 0;
    ssize_t length;
    while ((length = getline(&line, &allocated, stdin)) >= 0) {
        if (memchr(line, '\0', (size_t)length) || !parse_uint(line, &digit)
            || (digit == 0 && index != 0)
            || index == UINT64_MAX) {
            fprintf(stderr, "Invalid partial quotient at index %" PRIu64 "\n", index);
            free(line);
            return 1;
        }
        if (digit > maximum) {
            maximum = digit;
        }
        ring[index % capacity] = digit;
        if (index - first >= radius) {
            emit_window(ring, first, index, index - radius, threshold);
        }
        ++index;
    }
    free(line);
    if (ferror(stdin) || index == first) {
        fputs("Unreadable or empty input\n", stderr);
        return 1;
    }
    uint64_t pending = index - first < radius ? first : index - radius;
    for (; pending < index; ++pending) {
        emit_window(ring, first, index - 1, pending, threshold);
    }
    if (fflush(stdout) != 0) {
        return 1;
    }
    fprintf(stderr, "count=%" PRIu64 " maximum=%" PRIu64 "\n", index - first, maximum);
    return 0;
}
