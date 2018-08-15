#include "MyForm.h"
#include <array>

using namespace System;
using namespace System::Windows::Forms;
[STAThread]
int main(array<String^>^ argv)
{
	Application::EnableVisualStyles();
	Application::SetCompatibleTextRenderingDefault(false);
	Optimization::MyForm windowsForm;
	Application::Run(%windowsForm);
}