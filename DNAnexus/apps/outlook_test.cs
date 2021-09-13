using System;
using System.Windows.Forms;
using Outlook;

class Hello
{
    static void Main()
    {
        Console.WriteLine("Hello, World");
    }

    private void SetCurrentFolder()
    {
        string folderName = "Unresponsive";
        Outlook.Folder inBox = (Outlook.Folder)
            Application.ActiveExplorer().Session.GetDefaultFolder
            (Outlook.OlDefaultFolders.olFolderInbox);
        try
        {
            Application.ActiveExplorer().CurrentFolder = inBox.
                Folders[folderName];
            Application.ActiveExplorer().CurrentFolder.Display();
        }
        catch
        {
            MessageBox.Show("There is no folder named " + folderName +
                ".", "Find Folder Name");
        }
    }
}


