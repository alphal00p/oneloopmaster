set pagination off
set breakpoint pending on
break mi_process_init
run
watch theap_main.tld
disable 1
continue
bt 9
print theap_main.tld
print tld_main.thread_id
continue
bt 9
print theap_main.tld
print tld_main.thread_id
