CFLAGS=-Wall 
LIBS=-O2 -L/usr/X11R6/lib -lm -lpthread -lX11
EXE=NNSC

all:
	mkdir -p res
	g++ $(CFLAGS) NNSC.cpp -o $(EXE) $(LIBS)

test:
	./NNSC -i ./data/test_img.jpg -k 350 -m 0.075 -outm test_img_labelmap.png -outb test_img_border.png

test_list:
	./scripts/test_list.sh ./data/list_file.txt ./data/ 450 0.1
