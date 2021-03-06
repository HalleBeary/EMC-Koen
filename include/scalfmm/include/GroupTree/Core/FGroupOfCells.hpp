
// Keep in private GIT
#ifndef FGROUPOFCELLS_HPP
#define FGROUPOFCELLS_HPP

#include "Utils/FAssert.hpp"
#include "Utils/FAlignedMemory.hpp"
#include "Containers/FTreeCoordinate.hpp"
#include "../StarPUUtils/FStarPUDefaultAlign.hpp"

#include <list>
#include <functional>

/**
* @brief The FGroupOfCells class manages the cells in block allocation.
*/
template <class SymbolCellClass, class PoleCellClass, class LocalCellClass>
class FGroupOfCells {
    /** One header is allocated at the beginning of each block */
    struct alignas(FStarPUDefaultAlign::StructAlign) BlockHeader{
        MortonIndex startingIndex;
        MortonIndex endingIndex;
        int         numberOfCellsInBlock;
        int         idxGlobal;               ///< Global index of the group in the Octree
        bool        isMine;                  ///< true if I own the group
    };

protected:
    //< The size of the memoryBuffer
    size_t allocatedMemoryInByte;
    //< Pointer to a block memory
    unsigned char* memoryBuffer;

    //< Pointer to the header inside the block memory
    BlockHeader*    blockHeader;
    //< Pointer to the indexes table inside the block memory
    MortonIndex*    cellIndexes;
    //< Pointer to the cells inside the block memory
    SymbolCellClass*      blockCells;
#ifndef SCALFMM_SIMGRID_NODATA
    //< The multipole data
    PoleCellClass* cellMultipoles;
    //< The local data
    LocalCellClass* cellLocals;
#endif
    //< To kown if the object has to delete the memory
    bool deleteBuffer;

public:
    using multipole_t       = PoleCellClass;
    using local_expansion_t = LocalCellClass;
    using symbolic_data_t   = SymbolCellClass;

    FGroupOfCells()
        : allocatedMemoryInByte(0), memoryBuffer(nullptr),
          blockHeader(nullptr), cellIndexes(nullptr), blockCells(nullptr),
#ifndef SCALFMM_SIMGRID_NODATA
          cellMultipoles(nullptr), cellLocals(nullptr),
#endif
          deleteBuffer(false){
    }

    void reset(unsigned char* inBuffer, const size_t inAllocatedMemoryInByte,
               unsigned char* inCellMultipoles, unsigned char* inCellLocals){
        if(deleteBuffer){
            for(int idxCellPtr = 0 ; idxCellPtr < blockHeader->numberOfCellsInBlock ; ++idxCellPtr){
                (&blockCells[idxCellPtr])->~SymbolCellClass();
#ifndef SCALFMM_SIMGRID_NODATA
                (&cellMultipoles[idxCellPtr])->~PoleCellClass();
                (&cellLocals[idxCellPtr])->~LocalCellClass();
#endif
            }
            FAlignedMemory::DeallocBytes(memoryBuffer);
#ifndef SCALFMM_SIMGRID_NODATA
            FAlignedMemory::DeallocBytes(cellMultipoles);
            FAlignedMemory::DeallocBytes(cellLocals);
#endif
        }
        // Move the pointers to the correct position
        allocatedMemoryInByte = (inAllocatedMemoryInByte);
        memoryBuffer = (inBuffer);
        blockHeader         = reinterpret_cast<BlockHeader*>(inBuffer);
        inBuffer += sizeof(BlockHeader);
        cellIndexes   = reinterpret_cast<MortonIndex*>(inBuffer);
        inBuffer += (blockHeader->numberOfCellsInBlock*sizeof(MortonIndex));
        blockCells          = reinterpret_cast<SymbolCellClass*>(inBuffer);
        inBuffer += (sizeof(SymbolCellClass)*blockHeader->numberOfCellsInBlock);
        FAssertLF(size_t(inBuffer-memoryBuffer) == allocatedMemoryInByte);
#ifndef SCALFMM_SIMGRID_NODATA
        cellMultipoles = (PoleCellClass*)inCellMultipoles;
        cellLocals     = (LocalCellClass*)inCellLocals;
#endif
        deleteBuffer = (false);
    }

    /**
     * Init from a given buffer
     * @param inBuffer
     * @param inAllocatedMemoryInByte
     */
    FGroupOfCells(unsigned char* inBuffer, const size_t inAllocatedMemoryInByte,
                  unsigned char* inCellMultipoles, unsigned char* inCellLocals)
        : allocatedMemoryInByte(inAllocatedMemoryInByte), memoryBuffer(inBuffer),
          blockHeader(nullptr), cellIndexes(nullptr), blockCells(nullptr),
#ifndef SCALFMM_SIMGRID_NODATA
          cellMultipoles(nullptr), cellLocals(nullptr),
#endif
          deleteBuffer(false){
        // Move the pointers to the correct position
        blockHeader         = reinterpret_cast<BlockHeader*>(inBuffer);
        inBuffer += sizeof(BlockHeader);
        cellIndexes   = reinterpret_cast<MortonIndex*>(inBuffer);
        inBuffer += (blockHeader->numberOfCellsInBlock*sizeof(MortonIndex));
        blockCells          = reinterpret_cast<SymbolCellClass*>(inBuffer);
        inBuffer += (sizeof(SymbolCellClass)*blockHeader->numberOfCellsInBlock);
        FAssertLF(size_t(inBuffer-memoryBuffer) == allocatedMemoryInByte);
#ifndef SCALFMM_SIMGRID_NODATA
        cellMultipoles = (PoleCellClass*)inCellMultipoles;
        cellLocals     = (LocalCellClass*)inCellLocals;
#endif
    }

    /**
 * @brief FGroupOfCells
 * @param inStartingIndex first cell morton index
 * @param inEndingIndex last cell morton index + 1
 * @param inNumberOfCells total number of cells in the interval (should be <= inEndingIndex-inEndingIndex)
 */
    FGroupOfCells(const MortonIndex inStartingIndex, const MortonIndex inEndingIndex, const int inNumberOfCells)
        : allocatedMemoryInByte(0), memoryBuffer(nullptr), blockHeader(nullptr), cellIndexes(nullptr), blockCells(nullptr),
#ifndef SCALFMM_SIMGRID_NODATA
          cellMultipoles(nullptr), cellLocals(nullptr),
#endif
          deleteBuffer(true){
        FAssertLF((inEndingIndex-inStartingIndex) >= MortonIndex(inNumberOfCells));
        // Total number of bytes in the block
        const size_t memoryToAlloc = sizeof(BlockHeader) + (inNumberOfCells*sizeof(MortonIndex))
                + (inNumberOfCells*sizeof(SymbolCellClass));

        // Allocate
        FAssertLF(0 <= int(memoryToAlloc) && int(memoryToAlloc) < std::numeric_limits<int>::max());
        allocatedMemoryInByte = memoryToAlloc;
        memoryBuffer = (unsigned char*)FAlignedMemory::AllocateBytes<32>(memoryToAlloc);
        FAssertLF(memoryBuffer);
        memset(memoryBuffer, 0, memoryToAlloc);

        // Move the pointers to the correct position
        unsigned char* ptrBuff = memoryBuffer;
        blockHeader         = reinterpret_cast<BlockHeader*>(ptrBuff);
        ptrBuff += sizeof(BlockHeader);
        cellIndexes   = reinterpret_cast<MortonIndex*>(ptrBuff);
        ptrBuff += (inNumberOfCells*sizeof(MortonIndex));
        blockCells          = reinterpret_cast<SymbolCellClass*>(ptrBuff);
        ptrBuff += (sizeof(SymbolCellClass)*inNumberOfCells);
        FAssertLF(size_t(ptrBuff-memoryBuffer) == allocatedMemoryInByte);

        // Init header
        blockHeader->startingIndex         = inStartingIndex;
        blockHeader->endingIndex           = inEndingIndex;
        blockHeader->numberOfCellsInBlock  = inNumberOfCells;
        blockHeader->idxGlobal             = -1;
        blockHeader->isMine                = false;
#ifndef SCALFMM_SIMGRID_NODATA
        cellMultipoles = (PoleCellClass*)FAlignedMemory::AllocateBytes<32>(inNumberOfCells*sizeof(PoleCellClass));
        cellLocals     = (LocalCellClass*)FAlignedMemory::AllocateBytes<32>(inNumberOfCells*sizeof(LocalCellClass));
#endif
        for(int idxCell = 0 ; idxCell < inNumberOfCells ; ++idxCell){
#ifndef SCALFMM_SIMGRID_NODATA
            new (&cellMultipoles[idxCell]) PoleCellClass();
            new (&cellLocals[idxCell]) LocalCellClass();
#endif
            cellIndexes[idxCell] = -1;
        }
    }

    /** Call the destructor of cells and dealloc block memory */
    ~FGroupOfCells(){
        if(deleteBuffer){
            for(int idxCellPtr = 0 ; idxCellPtr < blockHeader->numberOfCellsInBlock ; ++idxCellPtr){
                (&blockCells[idxCellPtr])->~SymbolCellClass();
#ifndef SCALFMM_SIMGRID_NODATA
                (&cellMultipoles[idxCellPtr])->~PoleCellClass();
                (&cellLocals[idxCellPtr])->~LocalCellClass();
#endif
            }
            FAlignedMemory::DeallocBytes(memoryBuffer);
#ifndef SCALFMM_SIMGRID_NODATA
            FAlignedMemory::DeallocBytes(cellMultipoles);
            FAlignedMemory::DeallocBytes(cellLocals);
#endif
        }
    }

    /** Give access to the buffer to send the data */
    const unsigned char* getRawBuffer() const{
        return memoryBuffer;
    }

    /** The the size of the allocated buffer */
    size_t getBufferSizeInByte() const {
        return allocatedMemoryInByte;
    }

    /** Give access to the buffer to send the data */
    const PoleCellClass* getRawMultipoleBuffer() const{
#ifndef SCALFMM_SIMGRID_NODATA
        return cellMultipoles;
#else
        return nullptr;
#endif
    }

    /** Give access to the buffer to send the data */
    PoleCellClass* getRawMultipoleBuffer() {
#ifndef SCALFMM_SIMGRID_NODATA
        return cellMultipoles;
#else
        return nullptr;
#endif
    }

    /** The the size of the allocated buffer */
    size_t getMultipoleBufferSizeInByte() const {
        return sizeof(PoleCellClass)*blockHeader->numberOfCellsInBlock;
    }

    /** Give access to the buffer to send the data */
    LocalCellClass* getRawLocalBuffer(){
#ifndef SCALFMM_SIMGRID_NODATA
        return cellLocals;
#else
        return nullptr;
#endif
    }

    /** Give access to the buffer to send the data */
    const LocalCellClass* getRawLocalBuffer() const{
#ifndef SCALFMM_SIMGRID_NODATA
        return cellLocals;
#else
        return nullptr;
#endif
    }

    /** The the size of the allocated buffer */
    size_t getLocalBufferSizeInByte() const {
        return sizeof(LocalCellClass)*blockHeader->numberOfCellsInBlock;
    }

    /** To know if the object will delete the memory block */
    bool getDeleteMemory() const{
        return deleteBuffer;
    }

    /** The index of the fist cell (set from the constructor) */
    MortonIndex getStartingIndex() const {
        return blockHeader->startingIndex;
    }

    /** The index of the last cell + 1 (set from the constructor) */
    MortonIndex getEndingIndex() const {
        return blockHeader->endingIndex;
    }

    /** The number of cell (set from the constructor) */
    int getNumberOfCellsInBlock() const {
        return blockHeader->numberOfCellsInBlock;
    }

    /** The index of the current block in the global range **/
    int getIdxGlobal() const{
        return blockHeader->idxGlobal;
    }

    /* set the index of the current block */
    void setIdxGlobal(int idx){
        blockHeader->idxGlobal = idx;
    }

    // Return true if the current block is owned by the local tree
    bool isMine() const {
        return blockHeader->isMine;
    }

    // declare the block owned by the local tree
    void declare_mine(){
        blockHeader->isMine = true;
    }

    /** The size of the interval endingIndex-startingIndex (set from the constructor) */
    MortonIndex getSizeOfInterval() const {
        return MortonIndex(blockHeader->endingIndex-blockHeader->startingIndex);
    }

    /** Return true if inIndex should be located in the current block */
    bool isInside(const MortonIndex inIndex) const{
        return blockHeader->startingIndex <= inIndex && inIndex < blockHeader->endingIndex;
    }

    /** Return the idx in array of the cell */
    MortonIndex getCellMortonIndex(const int cellPos) const{
        FAssertLF(cellPos < blockHeader->numberOfCellsInBlock);
        return cellIndexes[cellPos];
    }

    /** Check if a cell exist (by binary search) and return it index */
    int getCellIndex(const MortonIndex cellIdx) const{
        int idxLeft = 0;
        int idxRight = blockHeader->numberOfCellsInBlock-1;
        while(idxLeft <= idxRight){
            const int idxMiddle = (idxLeft+idxRight)/2;
            if(cellIndexes[idxMiddle] == cellIdx){
                return idxMiddle;
            }
            if(cellIdx < cellIndexes[idxMiddle]){
                idxRight = idxMiddle-1;
            }
            else{
                idxLeft = idxMiddle+1;
            }
        }
        return -1;
    }

    /** Check if a cell exist (by binary search) and return it index */
    int getFistChildIdx(const MortonIndex parentIdx) const{
        int idxLeft = 0;
        int idxRight = blockHeader->numberOfCellsInBlock-1;
        while(idxLeft <= idxRight){
            int idxMiddle = (idxLeft+idxRight)/2;
            if((cellIndexes[idxMiddle]>>3) == parentIdx){
                while(0 < idxMiddle && (cellIndexes[idxMiddle-1]>>3) == parentIdx){
                    idxMiddle -= 1;
                }
                return idxMiddle;
            }
            if(parentIdx < (cellIndexes[idxMiddle]>>3)){
                idxRight = idxMiddle-1;
            }
            else{
                idxLeft = idxMiddle+1;
            }
        }
        return -1;
    }

    /** Return true if inIndex is located in the current block and is not empty */
    bool exists(const MortonIndex inIndex) const {
        return isInside(inIndex) && (getCellIndex(inIndex) != -1);
    }


    SymbolCellClass& getSymbolic(const int cellPos) {
        FAssertLF(cellPos < blockHeader->numberOfCellsInBlock);
        return blockCells[cellPos];
    }
    const SymbolCellClass& getSymbolic(const int cellPos) const {
        FAssertLF(cellPos < blockHeader->numberOfCellsInBlock);
        return blockCells[cellPos];
    }

    PoleCellClass& getMultipole(const int cellPos) {
        FAssertLF(cellPos < blockHeader->numberOfCellsInBlock);
        return cellMultipoles[cellPos];
    }
    const PoleCellClass& getMultipole(const int cellPos) const {
        FAssertLF(cellPos < blockHeader->numberOfCellsInBlock);
        return cellMultipoles[cellPos];
    }

    LocalCellClass& getLocalExpansion(const int cellPos) {
        FAssertLF(cellPos < blockHeader->numberOfCellsInBlock);
        return cellLocals[cellPos];
    }
    const LocalCellClass& getLocalExpansion(const int cellPos) const {
        FAssertLF(cellPos < blockHeader->numberOfCellsInBlock);
        return cellLocals[cellPos];
    }

    /** Allocate a new cell by calling its constructor */
    template<typename... CellConstructorParams>
    void newCell(const MortonIndex inIndex, const int id, CellConstructorParams... args){
        FAssertLF(isInside(inIndex));
        FAssertLF(!exists(inIndex));
        FAssertLF(id < blockHeader->numberOfCellsInBlock);
        new((void*)&blockCells[id]) SymbolCellClass(args...);
        cellIndexes[id] = inIndex;
    }

    /** Iterate on each allocated cells */
    template<class... Args>
    void forEachCell(std::function<void(SymbolCellClass*,PoleCellClass*,LocalCellClass*,Args...)> function, Args...args)
    {
        for(int idxCellPtr = 0 ; idxCellPtr < blockHeader->numberOfCellsInBlock ; ++idxCellPtr){
            function(&blockCells[idxCellPtr],
#ifndef SCALFMM_SIMGRID_NODATA
                     &cellMultipoles[idxCellPtr],
                     &cellLocals[idxCellPtr],
#else
                     nullptr,
                     nullptr,
#endif
                     args...);
        }
    }

    void forEachCell(std::function<void(SymbolCellClass*,PoleCellClass*,LocalCellClass*)> function)
    {
        for(int idxCellPtr = 0 ; idxCellPtr < blockHeader->numberOfCellsInBlock ; ++idxCellPtr){
            function(&blockCells[idxCellPtr],
#ifndef SCALFMM_SIMGRID_NODATA
                     &cellMultipoles[idxCellPtr],
                     &cellLocals[idxCellPtr]
#else
                     nullptr,
                     nullptr
#endif
                );
        }
    }

    /** Extract for implicit MPI */


    size_t extractGetSizeSymbUp(const std::vector<int>& cellsToExtract) const {
        return cellsToExtract.size() * (sizeof(SymbolCellClass) + sizeof(PoleCellClass));
    }


    void extractDataUp(const std::vector<int>& cellsToExtract,
                     unsigned char* outputBuffer, const size_t outputBufferSize) const {
        FAssertLF(outputBuffer || outputBufferSize == 0);
        size_t idxValue = 0;
        for(size_t idxEx = 0 ; idxEx < cellsToExtract.size() ; ++idxEx){
            const int idCell = cellsToExtract[idxEx];
            FAssertLF(idCell < blockHeader->numberOfCellsInBlock);
            memcpy(&outputBuffer[idxValue],
                   &blockCells[idCell],
                   sizeof(SymbolCellClass));
            idxValue += sizeof(SymbolCellClass);
            FAssertLF(idxValue <= outputBufferSize);
            memcpy(&outputBuffer[idxValue],
                   &cellMultipoles[idCell],
                   sizeof(PoleCellClass));
            idxValue += sizeof(PoleCellClass);
            FAssertLF(idxValue <= outputBufferSize);
        }
        FAssertLF(idxValue == outputBufferSize);
    }

    void restoreDataUp(const std::vector<int>& cellsToExtract,
                     const unsigned char* intputBuffer, const size_t inputBufferSize){
        size_t idxValue = 0;
        for(size_t idxEx = 0 ; idxEx < cellsToExtract.size() ; ++idxEx){
            const int idCell = cellsToExtract[idxEx];
            memcpy(&blockCells[idCell],
                   &intputBuffer[idxValue],
                   sizeof(SymbolCellClass));
            idxValue += sizeof(SymbolCellClass);
            memcpy(&cellMultipoles[idCell],
                   &intputBuffer[idxValue],
                   sizeof(PoleCellClass));
            idxValue += sizeof(PoleCellClass);
        }
        FAssertLF(idxValue == inputBufferSize);
    }

    size_t extractGetSizeSymbDown(const std::vector<int>& cellsToExtract) const {
        return cellsToExtract.size() * (sizeof(SymbolCellClass) + sizeof(LocalCellClass));
    }

    void extractDataDown(const std::vector<int>& cellsToExtract,
                     unsigned char* outputBuffer, const size_t outputBufferSize) const {
        size_t idxValue = 0;
        for(size_t idxEx = 0 ; idxEx < cellsToExtract.size() ; ++idxEx){
            const int idCell = cellsToExtract[idxEx];
            memcpy(&outputBuffer[idxValue],
                   &blockCells[idCell],
                   sizeof(SymbolCellClass));
            idxValue += sizeof(SymbolCellClass);
            memcpy(&outputBuffer[idxValue],
                   &cellLocals[idCell],
                   sizeof(PoleCellClass));
            idxValue += sizeof(PoleCellClass);
        }
        FAssertLF(idxValue == outputBufferSize);
    }

    void restoreDataDown(const std::vector<int>& cellsToExtract,
                     const unsigned char* intputBuffer, const size_t inputBufferSize){
        size_t idxValue = 0;
        for(size_t idxEx = 0 ; idxEx < cellsToExtract.size() ; ++idxEx){
            const int idCell = cellsToExtract[idxEx];
            memcpy(&blockCells[idCell],
                   &intputBuffer[idxValue],
                   sizeof(SymbolCellClass));
            idxValue += sizeof(SymbolCellClass);
            memcpy(&cellLocals[idCell],
                   &intputBuffer[idxValue],
                   sizeof(PoleCellClass));
            idxValue += sizeof(PoleCellClass);
        }
        FAssertLF(idxValue == inputBufferSize);
    }
};


#endif // FGROUPOFCELLS_HPP
