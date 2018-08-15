#pragma once
#include"DataManager.h"
#include"DotNetUtilities.h"
#include "Mvector.h"
#include "Mmatrix.h"
#include <stack>

namespace WindowsFormsApplication_cpp {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// WindowsForm 的摘要
	/// </summary>
	public ref class WindowsForm : public System::Windows::Forms::Form
	{
	public:
		WindowsForm(void)
		{
			InitializeComponent();
			dataManager = new DataManager();
			lastInputLength = -1;
		}

	protected:
		/// <summary>
		/// 清除任何使用中的資源。
		/// </summary>
		~WindowsForm()
		{
			if (components)
			{
				delete components;
			}
		}

	private: System::Windows::Forms::MenuStrip^  menuStrip2;
	private: System::Windows::Forms::ToolStripMenuItem^  FileToolStripMenuItem;

	private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel1;
	private: System::Windows::Forms::ToolStripMenuItem^  LoadVectorToolStripMenuItem;



	private: System::Windows::Forms::FlowLayoutPanel^  flowLayoutPanel1;
	private: System::Windows::Forms::FlowLayoutPanel^  flowLayoutPanel2;
	private: System::Windows::Forms::Label^  OutputLabel;
	private: System::Windows::Forms::TextBox^  Output;

	private: System::Windows::Forms::Label^  InputLabel;
	private: System::Windows::Forms::TextBox^  Input;
	private: System::Windows::Forms::Label^  VectorLabel;
	private: System::Windows::Forms::ListBox^  List;




	protected:















	protected:

	private:
		/// <summary>
			DataManager* dataManager;
			Mvector* mvtr;
			Mmatrix* mmtx;
			String^ userInput;
			int lastInputLength;
	private: System::Windows::Forms::OpenFileDialog^  openFileDialog1;
	private: System::Windows::Forms::OpenFileDialog^  openFileDialog2;
	private: System::Windows::Forms::ToolStripMenuItem^  LoadMatrixStripMenuItem1;

			 /// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// 此為設計工具支援所需的方法 - 請勿使用程式碼編輯器修改
		/// 這個方法的內容。
		/// </summary>
		void InitializeComponent(void)
		{
			this->menuStrip2 = (gcnew System::Windows::Forms::MenuStrip());
			this->FileToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->LoadVectorToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->LoadMatrixStripMenuItem1 = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->tableLayoutPanel1 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->flowLayoutPanel1 = (gcnew System::Windows::Forms::FlowLayoutPanel());
			this->InputLabel = (gcnew System::Windows::Forms::Label());
			this->Input = (gcnew System::Windows::Forms::TextBox());
			this->VectorLabel = (gcnew System::Windows::Forms::Label());
			this->List = (gcnew System::Windows::Forms::ListBox());
			this->flowLayoutPanel2 = (gcnew System::Windows::Forms::FlowLayoutPanel());
			this->OutputLabel = (gcnew System::Windows::Forms::Label());
			this->Output = (gcnew System::Windows::Forms::TextBox());
			this->openFileDialog1 = (gcnew System::Windows::Forms::OpenFileDialog());
			this->openFileDialog2 = (gcnew System::Windows::Forms::OpenFileDialog());
			this->menuStrip2->SuspendLayout();
			this->tableLayoutPanel1->SuspendLayout();
			this->flowLayoutPanel1->SuspendLayout();
			this->flowLayoutPanel2->SuspendLayout();
			this->SuspendLayout();
			// 
			// menuStrip2
			// 
			this->menuStrip2->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) { this->FileToolStripMenuItem });
			this->menuStrip2->Location = System::Drawing::Point(0, 0);
			this->menuStrip2->Name = L"menuStrip2";
			this->menuStrip2->Size = System::Drawing::Size(384, 24);
			this->menuStrip2->TabIndex = 1;
			this->menuStrip2->Text = L"menuStrip2";
			// 
			// FileToolStripMenuItem
			// 
			this->FileToolStripMenuItem->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(2) {
				this->LoadVectorToolStripMenuItem,
					this->LoadMatrixStripMenuItem1
			});
			this->FileToolStripMenuItem->Name = L"FileToolStripMenuItem";
			this->FileToolStripMenuItem->Size = System::Drawing::Size(38, 20);
			this->FileToolStripMenuItem->Text = L"File";
			// 
			// LoadVectorToolStripMenuItem
			// 
			this->LoadVectorToolStripMenuItem->Name = L"LoadVectorToolStripMenuItem";
			this->LoadVectorToolStripMenuItem->Size = System::Drawing::Size(152, 22);
			this->LoadVectorToolStripMenuItem->Text = L"Load Vector";
			this->LoadVectorToolStripMenuItem->Click += gcnew System::EventHandler(this, &WindowsForm::LoadVectorToolStripMenuItem_Click);
			// 
			// LoadMatrixStripMenuItem1
			// 
			this->LoadMatrixStripMenuItem1->Name = L"LoadMatrixStripMenuItem1";
			this->LoadMatrixStripMenuItem1->Size = System::Drawing::Size(152, 22);
			this->LoadMatrixStripMenuItem1->Text = L"Load Matrix";
			this->LoadMatrixStripMenuItem1->Click += gcnew System::EventHandler(this, &WindowsForm::toolStripMenuItem1_Click);
			// 
			// tableLayoutPanel1
			// 
			this->tableLayoutPanel1->ColumnCount = 2;
			this->tableLayoutPanel1->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent,
				50)));
			this->tableLayoutPanel1->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent,
				50)));
			this->tableLayoutPanel1->Controls->Add(this->flowLayoutPanel1, 1, 0);
			this->tableLayoutPanel1->Controls->Add(this->flowLayoutPanel2, 0, 0);
			this->tableLayoutPanel1->Dock = System::Windows::Forms::DockStyle::Fill;
			this->tableLayoutPanel1->Location = System::Drawing::Point(0, 24);
			this->tableLayoutPanel1->Name = L"tableLayoutPanel1";
			this->tableLayoutPanel1->RowCount = 1;
			this->tableLayoutPanel1->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 100)));
			this->tableLayoutPanel1->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Absolute, 20)));
			this->tableLayoutPanel1->Size = System::Drawing::Size(384, 338);
			this->tableLayoutPanel1->TabIndex = 2;
			this->tableLayoutPanel1->Paint += gcnew System::Windows::Forms::PaintEventHandler(this, &WindowsForm::tableLayoutPanel1_Paint);
			// 
			// flowLayoutPanel1
			// 
			this->flowLayoutPanel1->Controls->Add(this->InputLabel);
			this->flowLayoutPanel1->Controls->Add(this->Input);
			this->flowLayoutPanel1->Controls->Add(this->VectorLabel);
			this->flowLayoutPanel1->Controls->Add(this->List);
			this->flowLayoutPanel1->Location = System::Drawing::Point(195, 3);
			this->flowLayoutPanel1->Name = L"flowLayoutPanel1";
			this->flowLayoutPanel1->Size = System::Drawing::Size(186, 332);
			this->flowLayoutPanel1->TabIndex = 0;
			// 
			// InputLabel
			// 
			this->InputLabel->Anchor = System::Windows::Forms::AnchorStyles::None;
			this->InputLabel->AutoSize = true;
			this->InputLabel->Font = (gcnew System::Drawing::Font(L"Microsoft JhengHei", 9, static_cast<System::Drawing::FontStyle>((System::Drawing::FontStyle::Bold | System::Drawing::FontStyle::Italic)),
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(136)));
			this->InputLabel->Location = System::Drawing::Point(3, 0);
			this->InputLabel->Name = L"InputLabel";
			this->InputLabel->Size = System::Drawing::Size(41, 16);
			this->InputLabel->TabIndex = 0;
			this->InputLabel->Text = L"Input";
			// 
			// Input
			// 
			this->Input->Location = System::Drawing::Point(3, 19);
			this->Input->Multiline = true;
			this->Input->Name = L"Input";
			this->Input->Size = System::Drawing::Size(180, 158);
			this->Input->TabIndex = 1;
			this->Input->TextChanged += gcnew System::EventHandler(this, &WindowsForm::Input_TextChanged);
			// 
			// VectorLabel
			// 
			this->VectorLabel->Anchor = System::Windows::Forms::AnchorStyles::None;
			this->VectorLabel->AutoSize = true;
			this->VectorLabel->Font = (gcnew System::Drawing::Font(L"Microsoft JhengHei", 9, static_cast<System::Drawing::FontStyle>((System::Drawing::FontStyle::Bold | System::Drawing::FontStyle::Italic)),
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(136)));
			this->VectorLabel->Location = System::Drawing::Point(3, 180);
			this->VectorLabel->Name = L"VectorLabel";
			this->VectorLabel->Size = System::Drawing::Size(59, 16);
			this->VectorLabel->TabIndex = 2;
			this->VectorLabel->Text = L"Data List";
			// 
			// List
			// 
			this->List->FormattingEnabled = true;
			this->List->ItemHeight = 12;
			this->List->Location = System::Drawing::Point(3, 199);
			this->List->Name = L"List";
			this->List->Size = System::Drawing::Size(180, 124);
			this->List->TabIndex = 3;
			this->List->SelectedIndexChanged += gcnew System::EventHandler(this, &WindowsForm::List_SelectedIndexChanged);
			// 
			// flowLayoutPanel2
			// 
			this->flowLayoutPanel2->Controls->Add(this->OutputLabel);
			this->flowLayoutPanel2->Controls->Add(this->Output);
			this->flowLayoutPanel2->Location = System::Drawing::Point(3, 3);
			this->flowLayoutPanel2->Name = L"flowLayoutPanel2";
			this->flowLayoutPanel2->Size = System::Drawing::Size(186, 332);
			this->flowLayoutPanel2->TabIndex = 1;
			// 
			// OutputLabel
			// 
			this->OutputLabel->Anchor = System::Windows::Forms::AnchorStyles::None;
			this->OutputLabel->AutoSize = true;
			this->OutputLabel->Font = (gcnew System::Drawing::Font(L"Microsoft JhengHei", 9, static_cast<System::Drawing::FontStyle>((System::Drawing::FontStyle::Bold | System::Drawing::FontStyle::Italic)),
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(136)));
			this->OutputLabel->Location = System::Drawing::Point(3, 0);
			this->OutputLabel->Name = L"OutputLabel";
			this->OutputLabel->Size = System::Drawing::Size(52, 16);
			this->OutputLabel->TabIndex = 0;
			this->OutputLabel->Text = L"Output";
			// 
			// Output
			// 
			this->Output->Font = (gcnew System::Drawing::Font(L"PMingLiU", 9, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(136)));
			this->Output->Location = System::Drawing::Point(3, 19);
			this->Output->Multiline = true;
			this->Output->Name = L"Output";
			this->Output->ReadOnly = true;
			this->Output->Size = System::Drawing::Size(183, 313);
			this->Output->TabIndex = 1;
			this->Output->TextChanged += gcnew System::EventHandler(this, &WindowsForm::Output_TextChanged);
			// 
			// openFileDialog1
			// 
			this->openFileDialog1->FileName = L"openFileDialog1";
			this->openFileDialog1->FileOk += gcnew System::ComponentModel::CancelEventHandler(this, &WindowsForm::openFileDialog1_FileOk);
			// openFileDialog2
			// 
			this->openFileDialog2->FileName = L"openFileDialog2";
			this->openFileDialog2->FileOk += gcnew System::ComponentModel::CancelEventHandler(this, &WindowsForm::openFileDialog2_FileOk);
			// 
			// WindowsForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 12);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(384, 362);
			this->Controls->Add(this->tableLayoutPanel1);
			this->Controls->Add(this->menuStrip2);
			this->Name = L"WindowsForm";
			this->Text = L"VectorExample";
			this->Load += gcnew System::EventHandler(this, &WindowsForm::WindowsForm_Load);
			this->menuStrip2->ResumeLayout(false);
			this->menuStrip2->PerformLayout();
			this->tableLayoutPanel1->ResumeLayout(false);
			this->flowLayoutPanel1->ResumeLayout(false);
			this->flowLayoutPanel1->PerformLayout();
			this->flowLayoutPanel2->ResumeLayout(false);
			this->flowLayoutPanel2->PerformLayout();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion

private: System::Void WindowsForm_Load(System::Object^  sender, System::EventArgs^  e) {
}
private: System::Void LoadVectorToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e) 
{
	//開啟Dialog
	openFileDialog1->ShowDialog();
}
private: System::Void toolStripMenuItem1_Click(System::Object^  sender, System::EventArgs^  e) {
	openFileDialog2->ShowDialog();
}
private: System::Void Input_TextChanged(System::Object^  sender, System::EventArgs^  e)
{
	//當Input textbox中的輸入改變時，便會進入此函式
	//取得向量資料
	std::vector<Vector> vectors = dataManager->GetVectors();
	std::vector<Matrix> matrixs = dataManager->getMatrix();
	//判斷輸入自元為'\n'，並防止取到字串-1位置
	if (Input->Text->Length-1 >= 0 &&Input->Text[Input->Text->Length - 1] == '\n')
	{
		//將使用者輸入字串(在userInput中)，依空白作切割
		array<String^> ^userCommand = userInput->Split(' ');
		//字串比較，若指令為"print"的情況
		bool find = false;
		if (userCommand[0] == "print")
		{
			//定意輸出暫存
			String^ outputTemp = "";
			//透過for迴圈，從向量資料中找出對應變數
			for (unsigned int i = 0; i < vectors.size();i++)
			{
				//若變數名稱與指令變數名稱符合
				if (userCommand[1] == gcnew String(vectors[i].Name.c_str()))
				{
					//將輸出格式存入暫存
					find = true;
					outputTemp += "[";
					//將輸出資料存入暫存
					for (unsigned int j = 0; j<vectors[i].Data.size(); j++)
					{
						outputTemp += vectors[i].Data[j].ToString();
						if (j != vectors[i].Data.size() - 1)
							outputTemp += ",";
					}
					//將輸出格式存入暫存，並且換行
					outputTemp += "]" + Environment::NewLine;
					//輸出暫存資訊
					Output->Text += gcnew String(vectors[i].Name.c_str()) +" = "+ outputTemp;
					break;
				}
			}
			if (!find)
				Output->Text += "Can't find the vector!" + Environment::NewLine;
			else find = false;
		}
		//反之則判斷找不到指令
		else if (userCommand[0]=="V:") {
			int i = 1;
			//確認是否有該向量
			bool check = false;
			//存入四則運算中的運算符號
			std::stack<std::string> op_stack;
			//後敘式輸出
			std::vector<std::string> out_stack;
			//標註後序式位置是誰，1是向量資料 -1是運算符號 0是純量
			std::vector<int> index;
			//endl是字串結尾的意思
			while (userCommand[i] != "endl") {
				//如果輸入的指令在menu裡，表示為operation，存入stack裡
				if (mvtr->to_find(userCommand[i])) {
					std::string temp = "";
					MarshalString(userCommand[i], temp);
					op_stack.push(temp);
					check = true;
				}
				else if (userCommand[i] == ")") {
					//如果為)，表示在此符號前全丟入outputstack裡。
					for (int j = 0 ; !op_stack.empty(); j++)
						if (op_stack.top() == "(") {
							op_stack.pop();
							break;
						}
						else {
							out_stack.push_back(op_stack.top());
							index.push_back(-1);
							op_stack.pop();
						}
				}
				else {
					//表示去找向量，若有他的名字存入output stack裡
					std::string temp = "";
					MarshalString(userCommand[i], temp);
					for (unsigned int j = 0; j < vectors.size(); j++)
						if (temp == vectors[j].Name.c_str()) {
							out_stack.push_back(vectors[j].Name);
							index.push_back(1);
							check = true;
							break;
						}
					if (mvtr->isDigit(temp)) {
						//不是向量 確認是否為純量
						check = true;
						out_stack.push_back(temp);
						index.push_back(0);
					}
				}
				i++;
				if (!check) {
					Output->Text += "you input wrong things!!!" + Environment::NewLine;
					break;
				}
			}
			//確保運算符號堆疊全部清空
			if (!op_stack.empty())
				for (int j = 0;!op_stack.empty(); j++) {
					out_stack.push_back(op_stack.top());
					index.push_back(-1);
					op_stack.pop();
				}
			for (int j = 0; j < out_stack.size()&&check; j++) {
				Output->Text += gcnew String(out_stack.at(j).c_str())+" ";
			}
			//starting evulate mathematics operation
			std::string result = mvtr->mathematics_operation(vectors, out_stack,index);
			Output->Text += Environment::NewLine+gcnew String(result.c_str())+Environment::NewLine;
		}
		else if (userCommand[0] == "M:") {
			//使用陣列
			int i = 1;
			//標註後序式位置是誰，1是向量資料 -1是運算符號 0是純量
			std::vector<int> index;
			//儲存運算結果
			//Matrix result;
			//輸出
			std::vector<std::string> out;
			while (userCommand[i] != "endl") {
				std::string temp;
				MarshalString(userCommand[i],temp);
				out.push_back(temp);
				i++;
			}

			try {
				bool check = mmtx->Deal_outstack(out, matrixs, index);
			}
			catch (std::exception &e) {
				Output->Text += Environment::NewLine + gcnew String(e.what()) + Environment::NewLine;
			}
			

			//丟進Mmatrix_cpp裡處理 之後 out 為正確後序式輸出

			for (int j = 0; j < out.size(); j++) 
				Output->Text +=gcnew String(out[j].c_str())+" ";

			std::string result = mmtx->mathematics_operation(index, matrixs, out);
			out.clear();
			index.clear();
			Output->Text += Environment::NewLine + gcnew String(result.c_str()) + Environment::NewLine;
		}
		else if (userCommand[0]=="clr") {
			//清空符號
			Output->Text = "";
			Input->Text = "";
		}
		else
			Output->Text += "-Command not found-" + Environment::NewLine;
		userInput = "";
	}
	else
	{
		//將使用者輸入字串(在Text box中)，依'\n'作切割
		array<String^> ^userCommand = Input->Text->Split('\n');
		//並將最後一行，作為目前使用者輸入指令
		userInput = userCommand[userCommand->Length-1];
	}
}
private: System::Void openFileDialog1_FileOk(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e) 
{
	//在Dial	og按下OK便會進入此函式
	//字串轉制string^ to string
	std::string tempFileName;
	MarshalString(openFileDialog1->FileName, tempFileName);
	//將檔案路徑名稱傳入dataManager
	dataManager->SetFileName(tempFileName);
	//從讀取讀取向量資料
	if (dataManager->LoadVectorData())
	{
		//將VectorList中項目先做清除
		List->Items->Clear();	
		//取得所有向量資料
		std::vector<Vector> vectors = dataManager->GetVectors();

		for (unsigned int i = 0; i < vectors.size(); i++)
		{
			//將檔案名稱存入暫存
			std::string tempString = vectors[i].Name;
			//將輸出格式存入暫存
			tempString += " [";
			//將輸出資料存入暫存
			for (unsigned int j = 0; j<vectors[i].Data.size(); j++)
			{
				std::string scalarString = std::to_string(vectors[i].Data[j]);
				tempString += scalarString.substr(0, scalarString.size() - 5);
				if (j != vectors[i].Data.size() - 1)
					tempString += ",";
			}
			//將輸出格式存入暫存
			tempString += "]";
			//將項目加入VectorList中
			List->Items->Add(gcnew String(tempString.c_str()));
		}
		Output->Text += "-Vector datas have been loaded-" + Environment::NewLine;
	}
}
//當按下Load Matrix 時，就針對矩陣取值
private: System::Void openFileDialog2_FileOk(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e) {
	std::string tempName;
	MarshalString(openFileDialog2->FileName, tempName);
	dataManager->setName(tempName);

	if (dataManager->LoadData()) {
		List->Items->Clear();
		std::vector<Matrix> matrixs = dataManager->getMatrix();
		//指印出矩陣維度
		for (unsigned int i = 0; i < matrixs.size(); i++) {
			std::string tempString = matrixs[i].Name +" : row & column ";
			tempString += std::to_string(matrixs[i].row_len) +"  ";
			tempString += std::to_string(matrixs[i].col_len);
			List->Items->Add(gcnew String(tempString.c_str()));
		}
	}
	Output->Text += "-Matrix datas have been loaded-" + Environment::NewLine;
}
private: System::Void Output_TextChanged(System::Object^  sender, System::EventArgs^  e) {
}
private: System::Void tableLayoutPanel1_Paint(System::Object^  sender, System::Windows::Forms::PaintEventArgs^  e) {
}
private: System::Void List_SelectedIndexChanged(System::Object^  sender, System::EventArgs^  e) {
}
};
}
