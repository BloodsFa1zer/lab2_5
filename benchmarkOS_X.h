#include <mach/mach.h>
#include <mach/message.h>  // for mach_msg_type_number_t
#include <mach/kern_return.h>  // for kern_return_t
#include <mach/task_info.h>
#include <stdio.h>

int counter(void) {
    kern_return_t error;
    mach_msg_type_number_t outCount;
    mach_task_basic_info_data_t taskinfo;

    taskinfo.virtual_size = 0;
    outCount = MACH_TASK_BASIC_INFO_COUNT;
    error = task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&taskinfo, &outCount);
    if (error == KERN_SUCCESS) {
        // type is mach_vm_size_t
        printf("vsize = %llu\n", (unsigned long long)taskinfo.virtual_size);
        return 0;
    } else {
        printf("error %d\n", (int)error);
        return 1;
    }
}