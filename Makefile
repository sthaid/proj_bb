CC       = gcc
CFLAGS   = -g -O2 -Wall 
LDFLAGS  = -lm 

TARGET   = bb
SOURCES  = bb.c utils.c

OBJ := $(SOURCES:.c=.o)

$(TARGET): $(OBJ)
	$(CC) -o $@ $(OBJ) $(LDFLAGS)

clean:
	rm -f $(TARGET) $(OBJ)

