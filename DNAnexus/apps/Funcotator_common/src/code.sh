main() {
    java -jar ${DX_FS_ROOT}/dxExecutorWdl.jar workflow Inputs ${HOME} -traceLevel 1 -streamFiles PerFile 
}