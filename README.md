# DSP-HW1-HMM

## Instruction
1. Compile required programs
```
make
```
2. Train the HMM model
```
./train $iter $model_init_path$ $seq_path $output_model_path
```
- for example
```
./train 100 model_init.txt data/train_seq_01.txt model_01.txt

```
3. Test the HMM model
```
./test $models_list_path $seq_path $output_result_path
```
- for example
```
./test modellist.txt data/test_seq.txt result.txt
```