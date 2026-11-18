import os
import subprocess

# Trimmomatic 실행 함수
def run_trimmomatic(mode, input_files, output_files, adapter_file, threads=32):
    """
    Trimmomatic 작업을 실행합니다.

    Parameters:
        mode (str): "SE" 또는 "PE"
        input_files (list): 입력 FASTQ 파일 경로 (Single-End은 1개, Paired-End는 2개)
        output_files (list): 출력 파일 경로 (Single-End은 1개, Paired-End는 4개)
        adapter_file (str): Adapter 파일 경로
        threads (int): 사용할 스레드 수 (기본값: 32)
    """
    cmd = [
        "java", "-jar", "/usr/local/bin/trimmomatic.jar", mode,
        "-threads", str(threads),
        "-phred33",  # Quality 점수 형식
        *input_files,  # 입력 파일
        *output_files,  # 출력 파일
        f"ILLUMINACLIP:{adapter_file}:2:30:10",  # Adapter 설정
        "LEADING:20",  # 앞부분 Quality 필터
        "TRAILING:20",  # 뒷부분 Quality 필터
        "SLIDINGWINDOW:4:20",  # Sliding Window 설정
        "MINLEN:50"  # 최소 길이 필터
    ]

    # 실행 및 출력
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print("Completed:", " ".join(input_files))

# 작업 디렉토리 설정
input_dir = "./input_fastq"
output_dir = "./output_fastq"
adapter_file = "/usr/local/bin/TruSeq3-PE.fa"  # Adapter 파일 경로
os.makedirs(output_dir, exist_ok=True)

# 모든 샘플 파일 이름 자동 탐지 (1과 2로 짝을 이루는 Pair만 추출)
all_files = sorted(f for f in os.listdir(input_dir) if f.endswith(".fastq.gz"))
paired_samples = {}
for f in all_files:
    base_name = "_".join(f.split("_")[:-1])  # 마지막 숫자(1 또는 2)를 제외한 이름
    if base_name not in paired_samples:
        paired_samples[base_name] = []
    paired_samples[base_name].append(f)

# Trimmomatic 실행
for base_name, files in paired_samples.items():
    if len(files) != 2:
        print(f"Skipping {base_name}: Pair not complete ({files})")
        continue

    # 파일 정렬 및 경로 추가
    files.sort()  # '1'이 먼저, '2'가 나중에 오도록 정렬
    input_1 = os.path.join(input_dir, files[0])
    input_2 = os.path.join(input_dir, files[1])
    output_files = [
        os.path.join(output_dir, f"{base_name}_forward_paired.fastq.gz"),
        os.path.join(output_dir, f"{base_name}_forward_unpaired.fastq.gz"),
        os.path.join(output_dir, f"{base_name}_reverse_paired.fastq.gz"),
        os.path.join(output_dir, f"{base_name}_reverse_unpaired.fastq.gz"),
    ]

    # Trimmomatic 실행
    run_trimmomatic(
        mode="PE",
        input_files=[input_1, input_2],
        output_files=output_files,
        adapter_file=adapter_file,
        threads=32
    )
