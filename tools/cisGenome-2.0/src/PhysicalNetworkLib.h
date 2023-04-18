/* ----------------------------------------------------------------------- */
/*                                 Macro                                   */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */ 
/*                                 Struct                                  */
/* ----------------------------------------------------------------------- */ 

/* ----------------------------------------------------------------------- */ 
/*              NETWORKNODE: structure for network nodes                   */
/* ----------------------------------------------------------------------- */ 
struct NETWORKNODE
{
	/* node ID */
	int nID;
	int nNetIndex;

	/* node name */
	struct tagString *strName;

	/* special tag */
	int nTag;

	/* edges linked to the node */
	struct NODEEDGEINTERFACE *pEdgeList;

	/* double values */
	struct DOUBLEMATRIX *pDV;

	/* byte values */
	struct BYTEMATRIX *pBV;

	/* int values */
	struct INTMATRIX *pIV;

	/* string values */
	int nSVnum;
	struct tagString **vSV;

	/* pointer to next node */
	struct NETWORKNODE *pNext;
};

/* ----------------------------------------------------------------------- */ 
/*              NETWORKEDGE: structure for network edges                   */
/* ----------------------------------------------------------------------- */ 
struct NETWORKEDGE
{
	/* edge ID */
	int nEdgeID;

	/* node IDs */
	int nNodeID1;
	int nNodeID2;

	/* interaction direction: 0=no direction; 1=Node1-->Node2; 2=Node2-->Node1;
	     3=Node1<-->Node2 */
	int nDirection;

	/* interaction type: 0=unknown; 1=protein-protein; 2=protein-DNA; 3=DNA-protein;
		4=protein-RNA; 5=RNA-protein */
	int nInteractionType;
};

/* ----------------------------------------------------------------------- */ 
/*   NODEEDGEINTERFACE: structure for recording edges linked to a node     */
/* ----------------------------------------------------------------------- */ 
struct NODEEDGEINTERFACE
{
	/* edge ID */
	int nEdgeID;

	/* pointer to next interface element */
	struct NODEEDGEINTERFACE *pNext;
};

/* ----------------------------------------------------------------------- */ 
/*   NETWORKPHYSICAL: a physical network                                   */
/* ----------------------------------------------------------------------- */ 
struct NETWORKPHYSICAL
{
	/* Maximum Node number */
    int nMaxNodeNum;
	/* Real Node number */
    int nRealNodeNum;
	/* Head Node */
	struct NETWORKNODE *pHeadNode;
	/* Nodes */
	struct NETWORKNODE **vNodes;
	/* Edge number */
    int nEdgeNum;
	/* Edges */
	struct NETWORKEDGE **vEdges;
};

/* ----------------------------------------------------------------------- */ 
/*   NETWORKSHORTESTPATH: structure for recording network shortest paths   */
/* ----------------------------------------------------------------------- */ 
struct NETWORKSHORTESTPATH
{
	/* node number */
	int nNodeNum;

	/* node ID */
	struct INTMATRIX *pNodeID;

	/* neighboring node */
	struct INTMATRIX **vNeighborNode;

	/* shortest distance */
	struct DOUBLEMATRIX **vShortDist;
};

/* ----------------------------------------------------------------------- */
/*                              Declaration                                */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */ 
/* Network_LinkBINDEntrez_Main:                                            */
/* Link entrez gene id.                                                    */
/* ----------------------------------------------------------------------- */ 
int Network_LinkBINDEntrez_Main(char strBindFile[], char strEntrezFile[], 
						   char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/* Network_LinkBINDEntrez_GetGeneID:                                       */
/* Get entrez gene id.                                                     */
/* ----------------------------------------------------------------------- */ 
/* ----------------------------------------------------------------------- */ 
/* Network_LinkBINDEntrez_GetGeneID:                                       */
/* Get entrez gene id.                                                     */
/* ----------------------------------------------------------------------- */ 
int Network_LinkBINDEntrez_GetGeneID(int nMGI, char strMT[], char strMA[], 
			int nGeneNum, struct INTMATRIX *pGeneId,
			struct LONGMATRIX *pRNASortId, struct DOUBLEMATRIX *pRNASort,
			struct DOUBLEMATRIX *pRNAId, struct tagString **vRNAAcc, 
			struct LONGMATRIX *pProteinSortId, struct DOUBLEMATRIX *pProteinSort,
			struct DOUBLEMATRIX *pProteinId, struct tagString **vProteinAcc,
			struct LONGMATRIX *pGenomicSortId, struct DOUBLEMATRIX *pGenomicSort,
			struct DOUBLEMATRIX *pGenomicId, struct tagString **vGenomicAcc);

/* ----------------------------------------------------------------------- */ 
/* Network_GetSubNet_Main:                                                 */
/* Get subnetwork that contains given nodes.                               */
/* ----------------------------------------------------------------------- */ 
int Network_GetSubNet_Main(char strNetWorkFile[], char strNodeFile[], 
						   int nExportAnnotation, char strAnnotationFile[], 
						   char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/* Network_CreateOrthogloNet_Main:                                         */
/* Create an ortholog network.                                             */
/* ----------------------------------------------------------------------- */ 
int Network_CreateOrthogloNet_Main(char strNetWorkFile[], char strHomoloFile[], 
								   int nSrcSpecies, int nDestSpecies, 
								   char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/* Network_FindShortestPath_Main:                                          */
/* Find shortest path of a network.                                        */
/* ----------------------------------------------------------------------- */ 
int Network_FindShortestPath_Main(char strNetWorkFile[], char strAnnotationFile[],
								  char strSourceNodeFile[], char strDestNodeFile[], 
								  char strOutFile[], int nMaxIterNum);

/* ----------------------------------------------------------------------- */ 
/* Network_FindShortestPath:                                               */
/* Find shortest path.                                                     */
/* ----------------------------------------------------------------------- */ 
struct NETWORKSHORTESTPATH *Network_FindShortestPath(struct NETWORKPHYSICAL *pNet, 
					struct INTMATRIX *pDestNode, int nMaxIterNum);

/* ----------------------------------------------------------------------- */ 
/* Network_FindShortestPath_InitDist:                                      */
/* Init distance for shortest path finding.                                */
/* ----------------------------------------------------------------------- */ 
int Network_FindShortestPath_InitDist(struct NETWORKPHYSICAL *pNet, struct NETWORKSHORTESTPATH *pPath);

/* ----------------------------------------------------------------------- */ 
/* Network_FindShortestPath_ExportResults:                                 */
/* Export shortest path.                                                   */
/* ----------------------------------------------------------------------- */ 
int Network_FindShortestPath_ExportResults(char strFileName[], 
				struct NETWORKPHYSICAL *pNet, struct NETWORKSHORTESTPATH *pPath, 
				struct INTMATRIX *pSourceNode, struct INTMATRIX *pDestNode);

/* ----------------------------------------------------------------------- */ 
/* Network_FindShortestPath_V:                                             */
/* Find shortest path.                                                     */
/* ----------------------------------------------------------------------- */ 
struct NETWORKSHORTESTPATH *Network_FindShortestPath_V(struct NETWORKPHYSICAL *pNet, 
					struct INTMATRIX *pDestNode, int nRow, int nCol, double dMaxNodeValue,
					int nMaxIterNum);

/* ----------------------------------------------------------------------- */ 
/* Network_FindShortestPath_V_InitDist:                                    */
/* Init distance for shortest path finding.                                */
/* ----------------------------------------------------------------------- */ 
int Network_FindShortestPath_V_InitDist(struct NETWORKPHYSICAL *pNet, int nRow, int nCol,
				struct NETWORKSHORTESTPATH *pPath, double dMaxNodeValue);

/* ----------------------------------------------------------------------- */ 
/* Network_FindShortestPath_V_ExportResults:                               */
/* Export shortest path.                                                   */
/* ----------------------------------------------------------------------- */ 
int Network_FindShortestPath_V_ExportResults(char strFileName[], 
				struct NETWORKPHYSICAL *pNet, int nRow, int nCol, double dMaxNodeValue,
				struct NETWORKSHORTESTPATH *pPath, 
				struct INTMATRIX *pSourceNode, struct INTMATRIX *pDestNode);

/* ----------------------------------------------------------------------- */ 
/* Network_LoadFromSIF:                                                    */
/* Load a physical network from a SIF file.                                */
/* ----------------------------------------------------------------------- */ 
struct NETWORKPHYSICAL *Network_LoadFromSIF(char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/* Network_LoadAnnotation:                                                 */
/* Annotate network nodes from a file.                                     */
/* ----------------------------------------------------------------------- */
int Network_LoadAnnotation(struct NETWORKPHYSICAL *pNet, char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/* Network_InitNodeDV:                                                     */
/* Create DV vector for network nodes.                                     */
/* ----------------------------------------------------------------------- */
int Network_InitNodeDV(struct NETWORKPHYSICAL *pNet, int nHeight, int nWidth);

/* ----------------------------------------------------------------------- */ 
/* Network_InitNodeBV:                                                     */
/* Create BV vector for network nodes.                                     */
/* ----------------------------------------------------------------------- */
int Network_InitNodeBV(struct NETWORKPHYSICAL *pNet, int nHeight, int nWidth);

/* ----------------------------------------------------------------------- */ 
/* Network_InitNodeIV:                                                     */
/* Create IV vector for network nodes.                                     */
/* ----------------------------------------------------------------------- */
int Network_InitNodeIV(struct NETWORKPHYSICAL *pNet, int nHeight, int nWidth);


/* ----------------------------------------------------------------------- */ 
/* Network_InitNodeSV:                                                     */
/* Create SV vector for network nodes.                                     */
/* ----------------------------------------------------------------------- */
int Network_InitNodeSV(struct NETWORKPHYSICAL *pNet, int nSVnum);

/* ----------------------------------------------------------------------- */ 
/* Network_WriteToSIF:                                                     */
/* Write a physical network to a SIF file.                                 */
/* ----------------------------------------------------------------------- */ 
int Network_WriteToSIF(struct NETWORKPHYSICAL *pNet, char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/* Network_LoadValue_FromRefGene:                                          */
/* Load double values from a refgene database.                             */
/* ----------------------------------------------------------------------- */ 
int Network_LoadValue_FromRefGene(struct NETWORKPHYSICAL *pNet, int nRow, int nCol,
				struct tagRefGene **vRefGeneDatabase, int nRefGeneNum, 
				struct BYTEMATRIX *pHasValue, struct DOUBLEMATRIX *pSigValue,
				int nTakeAbsoluteValue, int nChooseMax);

/* ----------------------------------------------------------------------- */ 
/* Network_LoadValue_FromRefGene_ForShortestPath:                          */
/* Load double values from a refgene database.                             */
/* ----------------------------------------------------------------------- */ 
int Network_LoadValue_FromRefGene_ForShortestPath(struct NETWORKPHYSICAL *pNet, 
				int nRow, int nCol,
				struct tagRefGene **vRefGeneDatabase, int nRefGeneNum, 
				struct BYTEMATRIX *pHasValue, struct DOUBLEMATRIX *pSigValue,
				int nTakeAbsoluteValue, int nChooseMax, double dMaxNodeValue);

/* ----------------------------------------------------------------------- */ 
/* Network_DepthSearch_GetMaxScore:                                        */
/* Get maximun cumulative scores in the depth search.                      */
/* ----------------------------------------------------------------------- */ 
int Network_DepthSearch_GetMaxScore(struct NETWORKPHYSICAL *pNet, int nRow,
									int nCol, int nDepth);

/* ----------------------------------------------------------------------- */ 
/* Network_DepthSearch_GetMaxScore_Recursive:                              */
/* Perform recursive depth search.                                         */
/* ----------------------------------------------------------------------- */ 
int Network_DepthSearch_GetMaxScore_Recursive(struct NETWORKPHYSICAL *pNet,
		struct NETWORKNODE *pNode, int nRow, int nCol, int nDepth,
		int nLayer, struct tagString *pPastPath, double dPastScore, 
		struct DOUBLEMATRIX *pMaxScore, struct BYTEMATRIX *pMaxHas,
		struct tagString **vMaxPath);

/* ----------------------------------------------------------------------- */ 
/* NETWORKNODECREATE:                                                      */
/* Create an instance of NETWORKNODE.                                      */
/* ----------------------------------------------------------------------- */ 
struct NETWORKNODE *NETWORKNODECREATE();

/* ----------------------------------------------------------------------- */ 
/* NETWORKNODEDESTROY:                                                     */
/* Destroy an instance of NETWORKNODE.                                     */
/* ----------------------------------------------------------------------- */ 
void NETWORKNODEDESTROY(struct NETWORKNODE **ppNode);

/* ----------------------------------------------------------------------- */ 
/* NETWORKEDGECREATE:                                                      */
/* Create an instance of NETWORKEDGE.                                      */
/* ----------------------------------------------------------------------- */ 
struct NETWORKEDGE *NETWORKEDGECREATE();

/* ----------------------------------------------------------------------- */ 
/* NETWORKEDGEDESTROY:                                                     */
/* Destroy an instance of NETWORKEDGE.                                     */
/* ----------------------------------------------------------------------- */ 
void NETWORKEDGEDESTROY(struct NETWORKEDGE **ppEdge);

/* ----------------------------------------------------------------------- */ 
/* NODEEDGEINTERFACECREATE:                                                */
/* Create an instance of NODEEDGEINTERFACE.                                */
/* ----------------------------------------------------------------------- */ 
struct NODEEDGEINTERFACE *NODEEDGEINTERFACECREATE();

/* ----------------------------------------------------------------------- */ 
/* NODEEDGEINTERFACEDESTROY:                                               */
/* Destroy an instance of NODEEDGEINTERFACE.                               */
/* ----------------------------------------------------------------------- */ 
void NODEEDGEINTERFACEDESTROY(struct NODEEDGEINTERFACE **ppInterface);

/* ----------------------------------------------------------------------- */ 
/* NODEEDGEINTERFACELIST_CLEARALL:                                         */
/* Clear a list of NODEEDGEINTERFACE.                                      */
/* ----------------------------------------------------------------------- */ 
void NODEEDGEINTERFACELIST_CLEARALL(struct NODEEDGEINTERFACE **ppInterfaceList);

/* ----------------------------------------------------------------------- */ 
/* NETWORKPHYSICALCREATE:                                                  */
/* Create an instance of NETWORKPHYSICAL.                                  */
/* ----------------------------------------------------------------------- */ 
struct NETWORKPHYSICAL *NETWORKPHYSICALCREATE();

/* ----------------------------------------------------------------------- */ 
/* NETWORKPHYSICALDESTROY:                                                 */
/* Destroy an instance of NETWORKPHYSICAL.                                 */
/* ----------------------------------------------------------------------- */ 
void NETWORKPHYSICALDESTROY(struct NETWORKPHYSICAL **ppNetwork);

/* ----------------------------------------------------------------------- */ 
/* NETWORKSHORTESTPATHCREATE:                                              */
/* Create an object for performing shortest path search.                   */
/* ----------------------------------------------------------------------- */ 
struct NETWORKSHORTESTPATH *NETWORKSHORTESTPATHCREATE(int nNodeNum, int nTargetNum);

/* ----------------------------------------------------------------------- */ 
/* NETWORKSHORTESTPATHDESTROY:                                             */
/* Destroy a shortest path object.                                         */
/* ----------------------------------------------------------------------- */ 
void NETWORKSHORTESTPATHDESTROY(struct NETWORKSHORTESTPATH **ppNetPath);

/* ----------------------------------------------------------------------- */ 
/* NETWORKNODE_ADDEDGE:                                                    */
/* Add an edge to a network node.                                          */
/* ----------------------------------------------------------------------- */ 
int NETWORKNODE_ADDEDGE(struct NETWORKNODE **ppNode, int nNodeID, int nEdgeID);

/* ----------------------------------------------------------------------- */ 
/* NODEEDGEINTERFACELIST_ADDHEAD:                                          */
/* Add a node-edge interface object to a node-edge interface list.         */
/* ----------------------------------------------------------------------- */ 
int NODEEDGEINTERFACELIST_ADDHEAD(struct NODEEDGEINTERFACE **ppInterfaceList, 
		struct NODEEDGEINTERFACE *pInterface);