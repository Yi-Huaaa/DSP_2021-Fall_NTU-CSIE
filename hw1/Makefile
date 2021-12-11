.PHONY: all clean

ITER=100
TARGET=train test
ACC=acc

all: $(TARGET)

train: ./src/train.c
	gcc ./src/train.c -o train  -std=c99 -O2 -lm -Iinc -w

test: ./src/test.c
	gcc ./src/test.c -o test  -std=c99 -O2 -lm -Iinc -w

run:
	./train $(ITER) model_init.txt ./data/train_seq_01.txt model_01.txt
	./train $(ITER) model_init.txt ./data/train_seq_02.txt model_02.txt
	./train $(ITER) model_init.txt ./data/train_seq_03.txt model_03.txt
	./train $(ITER) model_init.txt ./data/train_seq_04.txt model_04.txt
	./train $(ITER) model_init.txt ./data/train_seq_05.txt model_05.txt
	./test modellist.txt ./data/test_seq.txt result.txt 
acc:
	gcc ./src/acc.c -o acc  -std=c99 -O2 -lm -Iinc -w
	./acc modellist.txt ./data/test_seq.txt result.txt ./data/test_lbl.txt

clean:
	$(RM) *.o
	$(RM) $(TARGET)
	$(RM) $(ACC)