#!/bin/bash
while getopts "n:c:o:m:t:" opt; do
    case ${opt} in
        c ) # run with CUDA Global Memory
            str_run_with="CUDA_Global"
            echo "Trial run with ${OPTARG}:"
            run_text="./harrisDetectorCuda -i ${OPTARG} -o ${str_run_with}_result_${OPTARG}"
            eval "make ; $run_text ; ./testDiffs referenceCuda.pgm ${str_run_with}_result_${OPTARG}"
            echo ""
            last_optarg=$OPTARG
            str_to_match="Device"

            if eval "$run_text ; ./testDiffs referenceCuda.pgm ${str_run_with}_result_${OPTARG}" | grep -q 'diffs'; then
                echo "Images do not match!"
                exit 1
            fi

        ;;
        m ) # run with CUDA Constant Memory
            str_run_with="CUDA_Constant"
            echo "Trial run with ${OPTARG}:"
            run_text="./harrisDetectorCudaConst -i ${OPTARG} -o ${str_run_with}_result_${OPTARG}"
            eval "make ; $run_text ; ./testDiffs referenceCudaConst.pgm ${str_run_with}_result_${OPTARG}"
            echo ""
            last_optarg=$OPTARG
            str_to_match="Device"

            if eval "$run_text ; ./testDiffs referenceCudaConst.pgm ${str_run_with}_result_${OPTARG}" | grep -q 'diffs'; then
                echo "Images do not match!"
                exit 1
            fi

        ;;
        t ) # run with CUDA Texture Memory
            str_run_with="CUDA_Texture"
            echo "Trial run with ${OPTARG}:"
            run_text="./harrisDetectorCudaTexture -i ${OPTARG} -o ${str_run_with}_result_${OPTARG}"
            eval "make ; $run_text ; ./testDiffs referenceCudaTexture.pgm ${str_run_with}_result_${OPTARG}"
            echo ""
            last_optarg=$OPTARG
            str_to_match="Device"

            if eval "$run_text ; ./testDiffs referenceCudaTexture.pgm ${str_run_with}_result_${OPTARG}" | grep -q 'diffs'; then
                echo "Images do not match!"
                exit 1
            fi

        ;;
        o ) # run with OpenMP
            str_run_with="OpenMP"
            echo "Trial run with ${OPTARG}:"
            run_text="./harrisDetectorOpenMP -i ${OPTARG} -o ${str_run_with}_result_${OPTARG}"
            eval "make ; $run_text ; ./testDiffs referenceOpenMP.pgm ${str_run_with}_result_${OPTARG}"
            echo ""
            last_optarg=$OPTARG
            str_to_match="OpenMP"

            if eval "$run_text ; ./testDiffs referenceOpenMP.pgm ${str_run_with}_result_${OPTARG}" | grep -q 'diffs'; then
                echo "Images do not match!"
                exit 1
            fi

        ;;
        n ) # run one of the above n times
            arr_name=(${last_optarg//./ })
            image_name=${arr_name[0]}
            results_file_name="${str_run_with}_time_results_${image_name}.txt"
            rm $results_file_name

            echo "Image: $last_optarg"
            echo "Starting $OPTARG runs"

            i=0
            t=0
            t_sum=0

            while [[ $i -ne $OPTARG ]]
            do
                t=$($run_text | awk -F': ' '/'${str_to_match}' processing time: /{print $2}' | awk '{split($0,a," "); print a[1]}')
                i=$(($i+1))
                t_sum=$(bc <<< "${t_sum}+${t}")
                echo "run $i: $t (ms)"
                printf "$t\n" >> $results_file_name
            done

            t_mean=$(bc -l <<<"${t_sum}/${OPTARG}")
            echo "Mean time: ${t_mean} (ms)"
        ;;
        \? ) # if none of the corret arguments are passed
            echo "Usage: cmd [-c] input_image_name [-t] input_image_name [-m] input_image_name [-o] input_image_name [-n] n_runs"
            exit 1
        ;;
        : ) # if an option does not get an argument that it needs
            echo "Invalid option: $OPTARG requires an argument"
            exit 1
        ;;
    esac
done